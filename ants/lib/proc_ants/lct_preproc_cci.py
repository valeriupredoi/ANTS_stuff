# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
"""
The following module is designed to support the generation of the CCI source
for pre-processing.  As such, it is non-general and will likely make
assumptions about the characteristics of the source.

"""
import ants
import ants.decomposition as decomp
import ants.fileformats.cover_mapping as cover_mapping
import iris
import numpy as np
from proc_ants import lct_preproc_igbp


def update_cci_metadata(cci_cube):
    """
    Fill and correct CCI metadata.

    Parameters
    ----------
    cci_cube : :class:`~iris.cube.Cube`
        CCI source cube.

    """
    # Correct grid definition.  Guess bounds is not quite what we want as the
    # edge cells bounds are not equally spaced - badly defined.
    ants.utils.cube.guess_horizontal_bounds(cci_cube)

    y_coord = cci_cube.coord(axis="y")
    bounds = y_coord.bounds.copy()
    bounds[0, 0], bounds[-1, -1] = 90, -90
    y_coord.bounds = bounds

    x_coord = cci_cube.coord(axis="x")
    bounds = x_coord.bounds.copy()
    bounds[0, 0], bounds[-1, -1] = -180, 180
    x_coord.bounds = bounds

    # Set basic flag value set
    set_flag_metadata(cci_cube)


def _igbp_cci_merge_operation(source, target):
    cci_snow_ice = -36
    cci_ocean = -45
    igbp_snow_ice = 15
    source.data = (source.data == igbp_snow_ice).astype("float16")
    result = source.regrid(
        target, ants.regrid.rectilinear.Linear(extrapolation_mode="linear")
    )
    data = np.empty(result.shape, dtype=target.dtype)
    data.fill(cci_ocean)
    # Set any snow/ice points according 50% coverage threshold.
    data[result.data > 0.5] = cci_snow_ice
    result.data = data
    return result


def merge_igbp(cci_source, igbp_source):
    """
    Merge the IGBP Antarctic with the CCI source.

    Parameters
    ----------
    cci_source : :class:`~iris.cube.Cube`
        CCI source cube.
    igbp_source : :class:`~iris.cube.Cube`
        IGBP source cube.

    Returns
    -------
    : :class:`~iris.cube.Cube`
        CCI source with Antartic data taken from the IGBP.

    """

    # Define the target we want to replace which is the Antarctica region
    # missing in the CCI source.
    min_lat = -60
    latc = iris.Constraint(latitude=lambda cell: cell < min_lat)
    target = cci_source.extract(latc)

    # Fix the igbp metadata
    igbp_source = lct_preproc_igbp.pre_process(igbp_source)

    # Extract to the region of the target
    extract_con = ants.ExtractConstraint(target)
    igbp_source = igbp_source.extract(extract_con)

    # Derive Antarctica data as derived from the igbp source.
    cci_missing = decomp.decompose(_igbp_cci_merge_operation, igbp_source, target)

    # merge igbp with the cci dataset.
    x, y = ants.utils.cube.horizontal_grid(cci_missing)
    minx, maxx = x.points.min(), x.points.max()
    miny, maxy = y.points.min(), y.points.max()
    slices = ants.utils.cube.get_slices(cci_source, [miny, maxy], [minx, maxx])[0]
    data = ants.utils.dask.deferred_data_update(
        cci_source.lazy_data(), cci_missing.lazy_data(), slices
    )
    return cci_source.copy(data)


def set_flag_metadata(cube):
    """
    Override CCI source flag_value and flag_meaning metadata.

    Parameters
    ----------
    cube : :class:`~iris.cube.Cube`
        Original CCI source before any pre-processing of flag values has
        occurred.

    """
    flag_values = [
        10,
        11,
        12,
        20,
        30,
        40,
        50,
        60,
        61,
        62,
        70,
        71,
        72,
        80,
        81,
        82,
        90,
        100,
        110,
        120,
        121,
        122,
        -126,
        -116,
        -106,
        -105,
        -104,
        -103,
        -96,
        -86,
        -76,
        -66,
        -56,
        -55,
        -54,
        -46,
        -45,
        -36,
        -47,
    ]
    flag_meanings = (
        "cropland_rainfed",
        "herbaceous_cover",
        "tree_or_shrub_cover",
        "cropland_irrigated",
        "mosaic_cropland",
        "mosaic_natural_vegetation",
        "tree_broadleaved_evergreen_closed_to_open",
        "tree_broadleaved_deciduous_closed_to_open",
        "tree_broadleaved_deciduous_closed",
        "tree_broadleaved_deciduous_open",
        "tree_needleleaved_evergreen_closed_to_open",
        "tree_needleleaved_evergreen_closed",
        "tree_needleleaved_evergreen_open",
        "tree_needleleaved_deciduous_closed_to_open",
        "tree_needleleaved_deciduous_closed",
        "tree_needleleaved_deciduous_open",
        "tree_mixed",
        "mosaic_tree_and_shrub",
        "mosaic_herbaceous",
        "shrubland",
        "shrubland_evergreen",
        "shrubland_dedicious",
        "grassland",
        "lichens_and_mosses",
        "sparse_vegetation",
        "sparse_tree",
        "sparse_shrub",
        "sparse_herbaceous",
        "tree_cover_flooded_fresh_or_brakish_water",
        "tree_cover_flooded_saline_water",
        "shrub_or_herbaceous_cover_flooded",
        "urban",
        "bare_areas",
        "consolidated_bare_areas",
        "unconsolidated_bare_areas",
        "water_bodies",
        "sea_ocean_water",
        "snow_and_ice",
        "resolved_lake",
    )

    cover_mapping.set_flag_arrays(cube, flag_values, flag_meanings)
