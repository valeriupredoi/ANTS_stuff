# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.coord_systems
import ants.fileformats._ugrid
import iris
import numpy as np
from ants import tests
from iris.coords import DimCoord


def simple_3d_time_varying(times=3):
    cube = geodetic(shape=(times, 5, 6))
    coord = iris.coords.DimCoord(
        np.arange(times, dtype="i8"), "time", units="hours since epoch"
    )
    cube.add_dim_coord(coord, 0)
    return cube


def simple_4d_with_hybrid_height():
    """
    air_temperature / (K)               (time: 3; model_level_number: 4; latitude: 5; longitude: 6)
         Dimension coordinates:
              time                           x                      -            -             -
              model_level_number             -                      x            -             -
              latitude                       -                      -            x             -
              longitude                      -                      -            -             x
         Auxiliary coordinates:
              level_height                   -                      x            -             -
              sigma                          -                      x            -             -
              surface_altitude               -                      -            x             x
         Derived coordinates:
              altitude                       -                      x            x             x

    """  # noqa
    cube = geodetic((3, 4, 5, 6))

    coord = iris.coords.DimCoord(
        np.arange(3, dtype="i8"), "time", units="hours since epoch"
    )
    cube.add_dim_coord(coord, 0)

    coord = iris.coords.DimCoord(
        np.arange(4, dtype="i8") + 10,
        "model_level_number",
        units="1",
        attributes={"positive": "up"},
    )
    cube.add_dim_coord(coord, 1)

    coord = iris.coords.AuxCoord(
        np.arange(4, dtype="i8") + 40, long_name="level_height", units="m"
    )
    cube.add_aux_coord(coord, 1)

    coord = iris.coords.AuxCoord(
        np.arange(4, dtype="i8") + 50, long_name="sigma", units="1"
    )
    cube.add_aux_coord(coord, 1)

    coord = iris.coords.AuxCoord(
        np.arange(5 * 6, dtype="i8").reshape(5, 6) + 100,
        long_name="surface_altitude",
        units="m",
    )
    cube.add_aux_coord(coord, [2, 3])

    factory = iris.aux_factory.HybridHeightFactory(
        delta=cube.coord("level_height"),
        sigma=cube.coord("sigma"),
        orography=cube.coord("surface_altitude"),
    )
    cube.add_aux_factory(factory)
    cube.rename("air_temperature")
    cube.units = "K"
    return cube


def simple_4d_with_hybrid_pressure():
    cube = geodetic((3, 4, 5, 6))

    coord = iris.coords.DimCoord(
        np.arange(3, dtype="i8"), "time", units="hours since epoch"
    )
    cube.add_dim_coord(coord, 0)

    coord = iris.coords.DimCoord(
        np.arange(4, dtype="i8") + 10, "model_level_number", units="1"
    )
    cube.add_dim_coord(coord, 1)

    coord = iris.coords.AuxCoord(
        np.arange(4, dtype="i8") + 40, long_name="level_pressure", units="Pa"
    )
    cube.add_aux_coord(coord, 1)

    coord = iris.coords.AuxCoord(
        np.arange(4, dtype="i8") + 50, long_name="sigma", units="1"
    )
    cube.add_aux_coord(coord, 1)

    coord = iris.coords.AuxCoord(
        np.arange(5 * 6, dtype="i8").reshape(5, 6) + 100,
        long_name="surface_air_pressure",
        units="Pa",
    )
    cube.add_aux_coord(coord, [2, 3])

    factory = iris.aux_factory.HybridPressureFactory(
        delta=cube.coord("level_pressure"),
        sigma=cube.coord("sigma"),
        surface_air_pressure=cube.coord("surface_air_pressure"),
    )
    cube.add_aux_factory(factory)
    cube.rename("air_temperature")
    cube.units = "K"
    return cube


def _add_metadata(cube, name, stash):
    if stash is not None:
        cube.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi(stash)
    if name is not None:
        cube.rename(name)


def _create_x_y_and_data(shape, xlim, ylim, with_bounds=True, data=None):
    shape = np.array(shape)
    assert shape.ndim == 1
    xdim = len(shape) - 1
    ydim = len(shape) - 2

    if data is None:
        data = np.arange(np.product(shape)).reshape(shape)
        data = data.astype("int32")

        if xlim[0] > xlim[1]:
            slices = [slice(None)] * shape.size
            slices[xdim] = slice(None, None, -1)
            data = data[tuple(slices)]
        if ylim[0] > ylim[1]:
            slices = [slice(None)] * shape.size
            slices[ydim] = slice(None, None, -1)
            data = data[tuple(slices)]

    x_bound = np.linspace(xlim[0], xlim[1], endpoint=True, num=shape[xdim] + 1)
    x_bounds = np.array([x_bound[:-1], x_bound[1:]]).T
    x_points = x_bounds.mean(axis=1)

    y_bound = np.linspace(ylim[0], ylim[1], endpoint=True, num=shape[ydim] + 1)
    y_bounds = np.array([y_bound[:-1], y_bound[1:]]).T
    y_points = y_bounds.mean(axis=1)

    if not with_bounds:
        x_bounds = y_bounds = None
    x = {
        "points": x_points,
        "bounds": x_bounds,
        "dim": xdim,
    }
    y = {
        "points": y_points,
        "bounds": y_bounds,
        "dim": ydim,
    }
    return x, y, data


def gen_curvilinear_cube(crs, shape=None, **kwargs):
    """
    Generate curvilinear lat-lon cube by translating 1D cube on the specified
    coordinate system.

    """

    def unflatten_bounds(flat_bnds):
        """
        Re-shape bounds from the flattened bounds array provided.

        Counterclockwise starting from the bottom left::
        3-2
        | |
        0-1

        """
        m, n = flat_bnds.shape
        shape = (m - 1, n - 1, 4)
        bounds = np.zeros(shape, flat_bnds.dtype)
        # Bounds used anticlockwise indexing.
        bounds[..., 0] = flat_bnds[:-1, :-1]
        bounds[..., 3] = flat_bnds[1:, :-1]
        bounds[..., 2] = flat_bnds[1:, 1:]
        bounds[..., 1] = flat_bnds[:-1, 1:]
        return bounds

    cube = gen_regular_cube(crs, shape=shape, **kwargs)
    x, y = ants.utils.cube.horizontal_grid(cube)

    cart_src_crs = crs.as_cartopy_crs()
    tgt_crs = iris.coord_systems.GeogCS(6371229.0)
    cart_tgt_crs = tgt_crs.as_cartopy_crs()

    # Derive x and y points in lat-lon crs
    xx_pnt, yy_pnt = np.meshgrid(x.points, y.points)
    xyz = cart_tgt_crs.transform_points(cart_src_crs, xx_pnt, yy_pnt)
    xx_pnt, yy_pnt = xyz[..., 0], xyz[..., 1]

    # Derive x and y bounds in lat-lon crs
    # - Create flat bounds arrays to provide meshgrid
    x_bounds = np.zeros(x.points.size + 1)
    y_bounds = np.zeros(y.points.size + 1)
    x_bounds[:-1] = x.bounds[:, 0]
    x_bounds[-1] = x.bounds[-1, 1]
    y_bounds[:-1] = y.bounds[:, 0]
    y_bounds[-1] = y.bounds[-1, 1]
    x_bounds, y_bounds = np.meshgrid(x_bounds, y_bounds)

    # - Reconstruct 2D bounds after transforming to the target crs.
    xyz = cart_tgt_crs.transform_points(cart_src_crs, x_bounds, y_bounds)
    x_bounds = unflatten_bounds(xyz[..., 0])
    y_bounds = unflatten_bounds(xyz[..., 1])

    # Remove original coords.
    cube.remove_coord(x)
    cube.remove_coord(y)

    # Add multi-dim coords to cube
    acrs = ants.coord_systems.CFCRS(tgt_crs)
    x_coord = iris.coords.AuxCoord(xx_pnt, bounds=x_bounds)
    ants.utils.coord.set_crs(x_coord, "x", acrs)
    y_coord = iris.coords.AuxCoord(yy_pnt, bounds=y_bounds)
    ants.utils.coord.set_crs(y_coord, "y", acrs)
    cube.add_aux_coord(x_coord, [0, 1])
    cube.add_aux_coord(y_coord, [0, 1])
    return cube


def gen_regular_cube(
    crs,
    shape=None,
    with_bounds=True,
    xlim=None,
    ylim=None,
    data=None,
    name=None,
    stash=None,
):
    cart_crs = crs.as_cartopy_crs()
    is_geodetic = "geodetic" in cart_crs.__class__.__name__.lower()
    if xlim is None:
        if is_geodetic:
            xlim = (-180, 180)
        else:
            xlim = cart_crs.x_limits
    if ylim is None:
        if is_geodetic:
            ylim = (-90, 90)
        else:
            ylim = cart_crs.y_limits
    if shape is None:
        shape = data.shape

    x, y, data = _create_x_y_and_data(
        shape, xlim, ylim, with_bounds=with_bounds, data=data
    )
    cube = iris.cube.Cube(data)

    acrs = ants.coord_systems.CFCRS(crs)
    x_coord = DimCoord(x["points"], bounds=x["bounds"])
    ants.utils.coord.set_crs(x_coord, "x", acrs)
    y_coord = DimCoord(y["points"], bounds=y["bounds"])
    ants.utils.coord.set_crs(y_coord, "y", acrs)

    cube.add_dim_coord(y_coord, y["dim"])
    cube.add_dim_coord(x_coord, x["dim"])

    ants.utils.cube.set_crs(cube, crs)
    ants.utils.cube.derive_circular_status(cube)
    _add_metadata(cube, name, stash)
    return cube


def geodetic(shape=None, north_pole_lat=90.0, north_pole_lon=0.0, **kwargs):
    """
    Helper function for generating a global lat-lon cube of specified shape.

    See also :func:`gen_regular_cube`

    """
    ellipsoid = iris.coord_systems.GeogCS(6371229.0)
    if north_pole_lat != 90.0 or north_pole_lon != 0.0:
        crs = iris.coord_systems.RotatedGeogCS(
            north_pole_lat, north_pole_lon, ellipsoid=ellipsoid
        )
    else:
        crs = ellipsoid
    return gen_regular_cube(crs, shape=shape, **kwargs)


def osgb(shape=None, **kwargs):
    """
    Return a cube covering the extent of the OSGB projection with the specified
    shape.

    See also :func:`gen_regular_cube`

    """
    crs = ants.coord_systems.OSGB.crs
    cube = gen_regular_cube(crs, shape=shape, **kwargs)
    return cube


def mesh_C4(load_data=False, constraint=iris.Constraint("sample_data")):
    """
    Return a C4 UGrid mesh cube.

    No provision for changing resolution.

    Current version of data file on disk has two data variables: 'sample_data'
    and 'additional_sample_data'.

    Parameters
    ----------
    load_data : boolean
        Whether to include a data payload.  Note that load_data=False implies
        the resulting cube will have a face_face_connectivity attributes,
        while load_data=True implies the cube will not have this attribute.
        This is a consequence of the current implementation (load_data=False
        means to load a standard cube generated by the LFRic mesh generator;
        load_data=True means to load a standard cube saved with ants.save,
        which currently does not include the face_face_connectivity
        attribute).

    constraint : :class:`iris.Constraint`
        Constraint for loading a particular data variable from the stock cube.
        Only valid if load_data is True.

    Returns
    -------

    : :class:`iris.cube.Cube` or :class:`iris.cube.CubeList`
        A single cube is returned when loading just the mesh
        (i.e. load_data=False).  Otherwise, a cubelist is returned.

    """
    if load_data:
        result = ants.fileformats._ugrid.load(
            tests.get_data_path(("stock", "data_C4.nc")), constraint
        )
    else:
        result = ants.fileformats._ugrid.load_mesh(
            tests.get_data_path(("stock", "mesh_C4.nc"))
        )

    return result
