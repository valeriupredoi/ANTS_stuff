# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import warnings
import weakref

import ants
import iris
import numpy as np
from ants.regrid.esmf import _LatLonExtractor

try:
    import ESMF

    _ESMF_IMPORT_ERROR = False
except Exception as _ESMF_IMPORT_ERROR:
    msg = " {}\nProceeding without capabilities provided by ESMPy (ESMF)."
    warnings.warn(msg.format(str(_ESMF_IMPORT_ERROR)))

# Dev note: ESMPy has the flexibility to use any values to represent masked
# data, and uses a second array to represent the mask.  Here, we're following
# the convention from the ESMPy tutorials of 1=>masked, 0=>unmasked (handily,
# since python booleans are really ints, this means the mask from a numpy
# masked array can be used directly).

# Dev note: ESMPy allows arbitrary dimension ordering for grids.  Ants usage
# is defined once, here, and any other references to the dimension ordering
# should use this defintion.  When referring to examples in ESMPy, note that
# they use "x" and "y" for these variable names:
[LONGITUDE_DIM, LATITUDE_DIM] = [0, 1]


class CachedRegridManager(object):
    """Caches ESMPy regrid objects to avoid re-computing weights.

    Only works for a regrid from a regular grid to an irregular mesh."""

    # Dev note: replace this with functools.lru_cache when we go to python
    # >=3.2
    def __init__(self):
        self._cache = weakref.WeakValueDictionary()

    def _create_regrid(self, srcfield, dstfield):
        # Dev note: To allow grid to grid or mesh to grid regrids, the
        # hardcoded indexing to the field masks needs to be updated.

        # mask[1] => elements since destination is a mesh:
        # http://www.earthsystemmodeling.org/esmf_releases/public/ESMF_7_1_0r/esmpy_doc/html/mesh.html#ESMF.api.mesh.Mesh.mask
        dst_mask_values = None
        dst_mask = dstfield.grid.mask[1]
        if np.any(dst_mask == 1):
            dst_mask_values = np.ndarray([1], dtype=np.int32)

        # mask[0] => centres since source is a grid:
        # http://www.earthsystemmodeling.org/esmf_releases/public/ESMF_7_1_0r/esmpy_doc/html/grid.html#ESMF.api.grid.Grid.mask
        # and
        # http://www.earthsystemmodeling.org/esmf_releases/public/ESMF_7_1_0r/esmpy_doc/html/StaggerLoc.html#ESMF.api.constants.StaggerLoc
        src_mask_values = None
        src_mask = srcfield.grid.mask[0]
        if np.any(src_mask == 1):
            src_mask_values = np.ndarray([1], dtype=np.int32)

        # This is the expensive call for which we want to cache the result.
        return ESMF.Regrid(
            srcfield,
            dstfield,
            src_mask_values=src_mask_values,
            dst_mask_values=dst_mask_values,
            regrid_method=ESMF.RegridMethod.CONSERVE,
            ignore_degenerate=False,
            unmapped_action=ESMF.UnmappedAction.IGNORE,
            norm_type=ESMF.api.constants.NormType.FRACAREA,
        )

    def _create_key(self, srcfield, dstfield):
        # Want to cache based on the underlying meshes/grids of
        # srcfield/dstfield (this avoids issues where data differences cause
        # the regridder to be uncached when it should be cached).  Notice that
        # ESMPy regridders are built using Field objects (Grid/Mesh + data)
        # rather than using grids directly, but it's only the underlying
        # Grid/Mesh that affects the weights.
        return str((srcfield.grid, dstfield.grid))

    def get_regrid(self, srcfield, dstfield):
        """
        Returns ESMPy Regrid object to regrid from srcfield to dstfield.

        If a cached ESMPy Regrid suitable for the srcfield and dstfield
        exists, that is the Regrid that's returned.  Otherwise, it's computed
        and added to the cache.

        Parameters
        ----------
        srcfield : :class:`ESMF.Field`
            Field containg the source Grid/Mesh and data to be regridded.
        dstfield : :class:`ESMF.Field`
            Field containg the target Grid/Mesh to which the srcfield will be
            regridded.

        Returns
        -------
        : :class:`ESMF.Regrid`
            Regrid object to call to perform the regridding.
        """
        if self._create_key(srcfield, dstfield) not in self._cache:
            regrid = self._create_regrid(srcfield, dstfield)
            self._cache[self._create_key(srcfield, dstfield)] = regrid
        else:
            regrid = self._cache[self._create_key(srcfield, dstfield)]
        return regrid

    def clear(self):
        """Empties the cache."""
        self._cache.clear()


class _UGridRegridder(object):
    """
    Regridding scheme for regridding from regular lat/lon grids to UGrid mesh.

    """

    # Leans heavily on ESMPy example:
    # http://www.earthsystemmodeling.org/esmf_releases/public/ESMF_7_1_0r/esmpy_doc/html/examples.html#regridding-from-grid-to-mesh
    # (and other examples from that page), as well as existing ants ESMPy
    # implementation.

    manager = CachedRegridManager()

    def __init__(self, source, target):
        self.srcgrid = _build_grid(source)
        self.dstgrid = _build_mesh(target)
        self.target = target

    def __call__(self, source):
        """
        Apply the interpolation weights to the source field.

        Parameters
        ----------
        source : :class:`~iris.cube.Cube`
           The cube to be regridded onto the target.  Needs to be on the same
           horizontal grid as the source used to initialise this class.

        Returns
        -------
        : class:`~iris.cube.Cube`
           Target cube with regridded data.

        """
        result = _create_result(source, self.target)
        self._reject_invalid_source(source)
        srcfield = ESMF.Field(
            self.srcgrid, name="srcfield", staggerloc=ESMF.StaggerLoc.CENTER
        )
        dstfield = ESMF.Field(
            self.dstgrid, name="dstfield", meshloc=ESMF.api.constants.MeshLoc.ELEMENT
        )

        if len(source.shape) == 2:
            srcfield.data[:, :] = source.data
            result.data = self._regrid_2d_field(srcfield, dstfield)
        else:
            # More complex case - need to break this into lat/lon layers and
            # then reassemble into right shape.  Easier to do at data
            # (i.e. numpy) level rather than cube level (all necessary cube
            # metadata already on result cube).
            def _regrid_layer(source_layer):
                srcfield.data[:, :] = source_layer
                return self._regrid_2d_field(srcfield, dstfield)

            result_2d = result.data.reshape(-1, result.shape[-1])
            source_2d = source.data.reshape(-1, source.shape[-2], source.shape[-1])

            for i in range(result_2d.shape[0]):
                result_2d[i] = _regrid_layer(source_2d[i])
            result.data = result_2d.reshape(result.data.shape)
        # If the source isn't global, or contains masked or NaN values, we
        # can't be certain that the output should not have NaNs:
        result_may_contain_nans = (
            np.any(np.isnan(source.data))
            or np.ma.is_masked(source.data)
            or not ants.utils.cube.is_global(source)
        )
        result_has_nans = np.any(np.isnan(result.data))
        if result_has_nans and not result_may_contain_nans:
            raise RuntimeError(f"Cube {result.name()} contains unexpected NaN values.")

        return result

    def _regrid_2d_field(self, srcfield, dstfield):
        dstfield.data[:] = np.NaN
        regrid = self.manager.get_regrid(srcfield, dstfield)
        regridded_field = regrid(
            srcfield, dstfield, zero_region=ESMF.api.constants.Region.SELECT
        )
        regridded_data = np.ma.masked_array(
            regridded_field.data.round(decimals=12),
            mask=regridded_field.grid.mask[1],
        )
        return regridded_data

    @staticmethod
    def _reject_invalid_source(source):
        """
        Rejects invalid source cube.

        These checks are in addition to those provided by ESMPy.

        Parameters
        ----------
        source: :class:`iris.cube.Cube`
            A regular lat/lon cube, rather than a UGrid cube.

        """
        # Check lat/lon are last two coordinates, in any order:
        lat_lon = (source.coord(axis="y"), source.coord(axis="x"))

        number_of_dimensions = len(source.shape)
        final_coord = source.coord(dimensions=(number_of_dimensions - 1))
        penultimate_coord = source.coord(dimensions=(number_of_dimensions - 2))

        if (final_coord not in lat_lon) or (penultimate_coord not in lat_lon):
            msg = (
                "Unsupported dimension order.  Expect latitude and "
                "longitude to be the final two coordinates, but "
                "found {} and {}".format(penultimate_coord.name(), final_coord.name())
            )
            raise RuntimeError(msg)


def _create_result(source, target):
    # Attributes that are needed to describe the mesh:
    mesh_attributes = [
        "mesh_topology",
        "face_face_connectivity",
        "face_node_connectivity",
        "nodes",
    ]

    if "Conventions" in source.attributes:
        source.attributes.pop("Conventions")

    contradictions = [
        attribute for attribute in source.attributes if attribute in mesh_attributes
    ]
    if contradictions:
        # May need to actually handle this in real world cases.  For now,
        # reject though:
        message = (
            "One or more attributes are defined on both regrid source "
            "and target: {}".format(contradictions)
        )
        raise ValueError(message)

    # Shape needs to be truncated since we're going from 2D lat/lon to 1D
    # lat/lon:
    shape = list(source.shape)[:-1]
    shape[-1] = target.shape[-1]
    result = _extrude_cube(target, shape, source.dtype)

    # Lat/lon handed in _extrude_cube; still need to reassemble additional
    # coordinates though:
    lat_lon_dim = len(result.shape) - 1
    for coord in source.aux_coords:
        dimensions = source.coord_dims(coord)
        try:
            if dimensions[0] < lat_lon_dim:
                result.add_aux_coord(coord, dimensions)
        except IndexError:
            """scalar coordinate"""
            result.add_aux_coord(coord)
    for coord in source.dim_coords:
        dimensions = source.coord_dims(coord)
        if dimensions[0] < lat_lon_dim:
            result.add_dim_coord(coord, dimensions)

    result.metadata = source.metadata
    # We know this isn't clobbering existing attributes due to contradictions
    # check above:
    for attribute in source.attributes:
        result.attributes[attribute] = source.attributes[attribute]
    mesh_attributes.append("Conventions")
    for attribute in mesh_attributes:
        if attribute in target.attributes:
            result.attributes[attribute] = target.attributes[attribute]

    return result


def _extrude_cube(cube, shape, dtype):
    """
    Extrudes a UGrid cube across extra dimensions.

    Cube must be a 1D representation of a plane - i.e. latitude and longitude
    as aux coords on a common dimension.

    Parameters
    ----------

    cube : :class:`iris.cube.Cube`
        Ugrid cube to be extruded, with latitude and longitude both as aux
        coords on the final dimension.
    shape : iterable
        Shape for the final cube.  This includes the existing 1D dimension as
        the last element of the iterable.

    Returns
    -------
    : :class:`iris.cube.Cube`
        Cube of the correct shape with horizontal coordinates from input cube,
        but no coordinates on other dimensions.  Data payload type is same as
        input cube, but data values are all 0.

    """
    if shape[-1] != cube.shape[-1]:
        raise ValueError(
            "Invalid shape.  Got {} expected {}".format(shape[-1], cube.shape[0])
        )

    data = np.zeros(shape, dtype=dtype)
    result = iris.cube.Cube(data)
    result.metadata = cube.metadata
    cube_horizontal_coord_dim = len(cube.shape) - 1
    result_horizontal_coord_dim = len(result.shape) - 1
    [
        result.add_aux_coord(coord, result_horizontal_coord_dim)
        for coord in cube.coords(dimensions=cube_horizontal_coord_dim)
    ]

    return result


def _build_grid(cube):
    # This is mostly lifted from the existing build_field method, but the grid
    # construction is decoupled from the field construction here (mostly
    # because the field construction is common whether the underlying
    # locations are a Grid or a Mesh).

    # Variable names here have mostly been left consistent with existing
    # regrid's build_field method to make integration easier in future.  We
    # could make different choices, but it'll make more work later.

    # Changes from build_field are:

    # 1. Addition of calculation of centre_lats/centre_lons (needed to regrid
    #    to Mesh ELEMENTS - ESMPy does not yet support regridding to Mesh
    #    NODES).
    # 2. The exception to keeping the build_field names are that lats/lons
    #    have been renamed as corner_lats/corner_lons to avoid potential
    #    confusion with centre_lats/centre_lons.
    # 3. Hard-coded staggering and coordsys for the values needed for grid to
    #    mesh.

    def _esmpy_zonal_mean_workaround(cube):
        # Workaround for zonal mean support (ESMPy cannot currently handle a
        # grid with cell edges that extend more than half way around the
        # sphere).
        #
        # concatenation will remove aux factories which is an accepted issue
        # for the circumstances where this workaround is required.
        x_coord = cube.coord(axis="x")
        if (
            ants.utils.cube.is_global(cube, x_axis_only=True)
            and x_coord.points.size == 1
        ):
            columns = 12
            column_width = 360.0 / columns
            lower_longitude = np.linspace(-180.0, 180.0 - column_width, columns)
            upper_longitude = lower_longitude + column_width
            longitude_points = lower_longitude + column_width / 2.0
            longitude_bounds = np.stack((lower_longitude, upper_longitude), axis=-1)
            cubes = iris.cube.CubeList()
            for column in range(columns):
                column_cube = cube.copy(cube.lazy_data())
                column_cube.coord(axis="x").points = longitude_points[column]
                column_cube.coord(axis="x").bounds = longitude_bounds[column]
                cubes.append(column_cube)
            cube = cubes.concatenate_cube()
            cube.coord(axis="x").circular = True
        return cube

    cube = _esmpy_zonal_mean_workaround(cube)

    # Need to construct an ESMPy grid.  This means we first need to get
    # the nodes set up:

    extractor = _LatLonExtractor(cube, "corner")
    corner_lats = extractor.get_latitude()
    corner_lons = extractor.get_longitude()

    cellDims = np.array([corner_lons.shape[0] - 1, corner_lats.shape[1] - 1])

    grid = ESMF.Grid(
        max_index=cellDims,
        staggerloc=[ESMF.StaggerLoc.CENTER],
        coord_sys=ESMF.api.constants.CoordSys.SPH_DEG,
    )

    # Now need to set up the face centres.  From:
    # http://www.earthsystemmodeling.org/esmf_releases/public/ESMF_7_1_0r/esmpy_doc/html/examples.html#create-a-2d-grid
    centre_lons = cube.coord(axis="X")
    centre_lats = cube.coord(axis="Y")
    len_lats = len(centre_lats.points)
    len_lons = len(centre_lons.points)

    gridXCenter = grid.get_coords(LONGITUDE_DIM)
    grid_lons = np.tile(centre_lons.points, (len_lats, 1))
    grid_lons = grid_lons.reshape(len_lats, len_lons)
    gridXCenter[...] = grid_lons
    gridYCenter = grid.get_coords(LATITUDE_DIM)
    grid_lats = centre_lats.points.repeat(len_lons)
    grid_lats = grid_lats.reshape(len_lats, len_lons)
    gridYCenter[...] = grid_lats

    # Allocate space for the vertices, ESMF wants the first coordinate to
    # be longitudes.
    grid.add_coords(staggerloc=ESMF.StaggerLoc.CORNER, coord_dim=0)
    # Get pointers to the ESMF coordinates.
    lonPoint = grid.get_coords(coord_dim=0, staggerloc=ESMF.StaggerLoc.CORNER)
    latPoint = grid.get_coords(coord_dim=1, staggerloc=ESMF.StaggerLoc.CORNER)

    ibeg0 = grid.lower_bounds[ESMF.StaggerLoc.CORNER][0]
    iend0 = grid.upper_bounds[ESMF.StaggerLoc.CORNER][0]
    ibeg1 = grid.lower_bounds[ESMF.StaggerLoc.CORNER][1]
    iend1 = grid.upper_bounds[ESMF.StaggerLoc.CORNER][1]
    lonPoint[...] = corner_lons[ibeg0:iend0, ibeg1:iend1]
    latPoint[...] = corner_lats[ibeg0:iend0, ibeg1:iend1]

    # Now add the mask to the CENTER (default) staggerloc:
    mask = grid.add_item(ESMF.GridItem.MASK)
    mask[:] = 0
    if np.ma.is_masked(cube.data):
        latitude_name = cube.coord(axis="y").name()
        longitude_name = cube.coord(axis="x").name()
        layer = next(cube.slices([latitude_name, longitude_name]))
        mask[layer.data.mask] = 1

    return grid


def _build_mesh(cube):
    """
    Builds ESMPy mesh constructed from the cube.

    Assumes cube represents UGRID data.

    """
    num_elem = cube.shape[-1]

    # TODO: mapping_coord should be derived from the cube
    face_node_mapping = cube.attributes["face_node_connectivity"].data

    # First create a list of nodes expressed as (node index, lon, lat),
    # and use a set to eliminate duplicates.  Note the swapping of the
    # lat/lon coords:
    node_coord = set(
        [
            (index, lon, lat)
            for (index, lat, lon) in zip(
                face_node_mapping.flatten(),
                cube.coord("latitude").bounds.flatten(),
                cube.coord("longitude").bounds.flatten(),
            )
        ]
    )

    # Then sort and convert to numpy array, while dropping out the index
    # coordinate.
    node_coord = np.array(list(sorted(node_coord)))[:, 1:]

    num_node = node_coord.shape[0]

    # Index for nodes:
    nodeID = np.arange(num_node)

    # All nodes on same processor, i.e. 0: we don't want to use ESMPy
    # parallelism
    nodeOwner = np.zeros(num_node).flatten()

    # Index for elements:
    elemID = np.arange(num_elem)

    # All elements are quadrilateral:
    elemType = np.ones_like(elemID) * ESMF.MeshElemType.QUAD
    elemType = elemType.flatten()

    # Associating the correct nodes to the element.  Note the -1 since the
    # mapping coordinate is 1-indexed:
    #
    # TODO: UGrid spec does allow mapping to be either zero or one based,
    # with the start_index attribute set:
    # http://ugrid-conventions.github.io/ugrid-conventions/#zero-or-one-based-indexing
    # We should check this attribute and handle it accordingly.
    elemConn = face_node_mapping.flatten() - 1

    # And the element locations:
    elemCoord = np.array(
        [
            cube.coord("longitude").points.flatten(),
            cube.coord("latitude").points.flatten(),
        ]
    ).T

    # Finally construct the mesh from the node and element info:
    mesh = ESMF.Mesh(
        parametric_dim=2, spatial_dim=2, coord_sys=ESMF.api.constants.CoordSys.SPH_DEG
    )
    mesh.add_nodes(num_node, nodeID, np.array(node_coord).flatten(), nodeOwner)

    # Mask needs to be added with the elements (contrast with build_grid,
    # where mask is added after the fact)
    element_mask = np.zeros(num_elem)
    if np.ma.is_masked(cube.data):
        if len(cube.shape) == 1:
            # Cube is a single lat/lon layer
            element_mask = cube.data.mask
        else:
            # Only need single layer for the mask
            index = np.zeros(cube.shape, dtype="bool")
            index[0] = True
            element_mask[cube.data.mask[index]] = 1

    mesh.add_elements(
        num_elem,
        elemID,
        elemType,
        elemConn,
        element_coords=elemCoord,
        element_mask=element_mask,
    )
    return mesh


class _UGrid(object):
    """
    ESMPy regridder for mesh data.

    Currently, only supports lat/lon grid to UGrid mesh.

    """

    def regridder(self, source_grid, target_mesh):
        """
        Assumes source data is on cell centres, and regrids to mesh elements,
        using the ESMPy conservative regridding method.  These constraints are
        flexible and can be relaxed if there's a use case.

        Parameters
        ----------
        source_grid : :class:`~iris.cube.Cube`
            Regular lat/lon cube used to define the source grid.
        target_mesh : :class:`~iris.cube.Cube`
            UGrid mesh cube used to define the target mesh.

        Returns
        -------
        : :class:`ESMFRegridder`
           Callable with the interface `callable(cube)`

           where `cube` is a cube with the same grid as `source_grid`
           that is to be regridded to the `target_mesh`.

        """
        return _UGridRegridder(source_grid, target_mesh)
