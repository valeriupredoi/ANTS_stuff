# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import itertools
import logging
import warnings

import ants
import cartopy.io.shapereader as shpreader
import numpy as np
import shapely
from shapely.vectorized import contains

_LOGGER = logging.getLogger(__name__)


def get_lake_geoms():
    def get_geom(res, name, category):
        filename = shpreader.natural_earth(resolution=res, name=name, category=category)
        lakes = shpreader.Reader(filename)
        records = lakes.records()
        return records

    geom1 = get_geom("10m", "lakes", "physical")
    geom2 = get_geom("50m", "lakes", "physical")
    geom3 = get_geom("110m", "lakes", "physical")
    return itertools.chain(geom1, geom2, geom3)


def get_lake_geom(name):
    records = get_lake_geoms()
    res = None
    _LOGGER.info(f"Attempting to find geometry for {name}")
    for record in records:
        if name in record.attributes["name"].lower():
            res = record
            break
    if not res:
        warnings.warn(f"Unable to find Natural Earth shapefile for {name}")
    return res


def extract_region(cube, geom):
    """
    Given a cube and a geometry, extract the cube sub-region corresponding to
    that geometry.
    """
    # Get generous bounding box for geometry (lake)
    #   - We buffer this geometry by a percentage of the region extent.
    # containment_geom is the geometry we use to constrain the floodfill.
    # extraction_geom is the geometry we use to extract the cube sub region.
    minx, miny, maxx, maxy = geom.bounds
    x = np.array([minx, minx, maxx, maxx])
    y = np.array([miny, maxy, miny, maxy])
    distances = np.sqrt((x[1:] - x[:-1]) ** 2 + (y[1:] - y[:-1]) ** 2)
    distance = distances.max()
    containment_geom = geom.buffer(distance * 0.02)
    extraction_geom = geom.buffer(distance * 0.25)

    # Slice out region corresponding to geometry (lake)
    minx, miny, maxx, maxy = extraction_geom.bounds
    slices = ants.utils.cube.get_slices(cube, [miny, maxy], [minx, maxx])
    if len(slices) != 1:
        msg = "Expecting 1 slice, got {} slices".format(len(slices))
        raise RuntimeError(msg)
    slices = slices[0]
    region_cube = cube.copy(cube.lazy_data())[slices]
    region_cube.data = np.ma.array(region_cube.data.astype("int8"))
    return region_cube, extraction_geom, containment_geom, slices


def fetch_seed_index(cube, seed):
    """
    Given seed value in lat-lon, return the corresponding index within the
    cube provided.
    """
    x, y = ants.utils.cube.horizontal_grid(cube, dim_coords=True)
    if (
        seed[1] < x.points.min()
        or seed[1] > x.points.max()
        or seed[0] < y.points.min()
        or seed[0] > y.points.max()
    ):
        msg = (
            "Seed value x,y:{} is not contained within the extent of the "
            "extracted domain: xlim [{}, {}], ylim [{}, {}]"
        )
        raise ValueError(
            msg.format(
                seed,
                x.points.min(),
                x.points.max(),
                y.points.min(),
                y.points.max(),
            )
        )
        return cube
    xd = abs(x.points - seed[1]).argmin()
    yd = abs(y.points - seed[0]).argmin()
    return xd, yd


def fill_lakes(
    cube,
    lake_name,
    seed,
    src_flag_meanings,
    tgt_flag_meaning,
    constrain=False,
    use_lake_geom=True,
):
    """
    Fill source type with the target type specified

    Parameters
    ----------
    cube : :class:`iris.cube.Cube`
        Source cube to which we want to have the specified source filled.
    lake_name : str
        Lake name used to constrain the cube to the relevant area.
    seed : tuple
        Floodfill seed points for both lat and lon in units of degrees of form
        tuple(y, x).
    src_flag_meanings : str
        Source type name(s) lookup for the provided source cube.
    tgt_flag_meaning : str
        Target type name lookup for the provided source cube.
    constrain : bool
        Determine whether the floodfill should be constrained by the lake
        geometry.
    use_lake_geom : bool
        Use Natural Earth gometry to extract region.  Note that 'constrain' will
        determine whether the floodfill is subsequently constrained by this
        geometry.

    Returns
    -------
    : :class:`iris.cube.Cube`
        Source classification data.

    """
    if isinstance(src_flag_meanings, str):
        src_flag_meanings = [src_flag_meanings]

    flag_values, flag_meanings = ants.fileformats.cover_mapping.get_flag_arrays(cube)
    # Fetch flag values corresponding to the src flag strings provided.
    inds = [
        flag_meanings.tolist().index(src_flag_meaning.lower())
        for src_flag_meaning in src_flag_meanings
    ]
    src_values = [int(flag_values[ind]) for ind in inds]

    # Fetch flag value corresponding to the tgt flag string provided.
    ind = flag_meanings.tolist().index(tgt_flag_meaning.lower())
    tgt_value = int(flag_values[ind])

    # Fetch the geometry for the lake name specified.
    geom = None
    if use_lake_geom:
        geom = get_lake_geom(lake_name)
    if geom is None:
        # No Natural earth geometry found so lets pull out a region.
        width = 0.8
        minlon, maxlon = seed[1] - width, seed[1] + width
        minlat, maxlat = seed[0] - width, seed[0] + width
        geom = shapely.geometry.Polygon(
            (
                (minlon, minlat),
                (minlon, maxlat),
                (maxlon, maxlat),
                (maxlon, minlat),
                (minlon, minlat),
            )
        )
    else:
        geom = geom.geometry

    lake_cube, extraction_geom, containment_geom, slices = extract_region(cube, geom)
    lake_data = lake_cube.data

    if 1 in lake_cube.shape:
        warnings.warn(
            f"Lake {lake_name} too small to be filled at this resolution "
            f"(seed: {seed}). Skipping."
        )
        return cube
    xd, yd = fetch_seed_index(lake_cube, seed)

    # Perform floodfill operation
    seed_value = lake_data[yd, xd]
    if seed_value in src_values:
        if constrain:
            # Mask out region not overlapping lake geometry to constrain floodfill.
            x, y = ants.utils.cube.horizontal_grid(lake_cube, dim_coords=True)
            xx, yy = np.meshgrid(x.points, y.points)
            mask_not_river = ~contains(containment_geom, xx, yy) * ~np.ma.getmaskarray(
                lake_data
            )
            lake_data.mask = mask_not_river
        ants.analysis.floodfill(lake_data, (yd, xd), tgt_value)
        if constrain:
            # Remove the constraint (mask).
            lake_data.mask[mask_not_river] = False

        data = ants.utils.dask.deferred_data_update(cube.lazy_data(), lake_data, slices)
        # Creating a copy is neccessary as the data is private and we cannot
        # subset a view of the cube.
        # https://github.com/SciTools/iris/pull/1992
        cube = cube.copy(data)
    elif seed_value == tgt_value:
        msg = """Data already has target type "{}" at seed {} for """
        msg += "lake {}"
        warnings.warn(msg.format(tgt_flag_meaning, seed, lake_name))
    else:
        msg = """Data doesn't have "{}" type at seed {} for lake {}"""
        warnings.warn(msg.format(src_flag_meanings, seed, lake_name))

    return cube
