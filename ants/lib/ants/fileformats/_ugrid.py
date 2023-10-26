# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import warnings
from collections import namedtuple

import ants
import dask.array as da
import iris
import numpy as np
from netCDF4 import Dataset

# Define the cf_roles that are special cases on both load and save
UGRID_ROLES = ["mesh_topology", "face_face_connectivity", "face_node_connectivity"]


def load(filename, data_constraints=None):
    """
    Loads a UGRID file including data variables.

    Expectation is that the file not only conforms to the UGRID specification,
    but also conforms to XIOS compatible output.  In particular, this means
    that the variables in the UGRID file are expected to have long names
    which are requirements over and above the UGRID specification.

    Parameters
    ----------

    filename : str
        Name of the UGRID file to load.

    data_constraint : list of :class:`iris.Constraint` (optional)
        Iris constraint to identify which data variables to load.  If omitted,
        all data variables will be loaded.

    Returns
    -------

    : :class:`iris.cube.CubeList`
        Each cube has a data payload, and 1D coordinates for each of latitude
        and longitude, with the nodes from the original mesh translated into
        4-element bounds.  The mesh topology, the face to node, and the face
        to face mapping are attached as attributes.

    """
    cubes = iris.load(filename)
    if data_constraints:
        potential_data_cubes = cubes.extract_strict(data_constraints)
        # extract_strict returns single cube for a single constraint, or
        # multiple cubes for no or multiple constraints:
        potential_data_cubes = ants.utils.cube.as_cubelist(potential_data_cubes)
    else:
        potential_data_cubes = cubes
    data_variables = [
        c.name()
        for c in potential_data_cubes
        if "mesh" in c.attributes and "location" in c.attributes
    ]
    result = iris.cube.CubeList()
    for data_variable in data_variables:
        ugrid = _UGridCubes(cubes, data_variable)
        result.append(ugrid._construct_mesh())

    # Note that iris.load with a single constraint still returns a cubelist,
    # (unlike extract_strict) so the return from this load function is always
    # a cubelist for consistency:
    return result


def load_mesh(filename):
    """
    Loads a UGRID mesh from a file.

    Only the mesh is loaded - any data variables are ignored.

    Expectation is that the file not only conforms to the UGRID specification,
    but also conforms to XIOS compatible output.  In particular, this means
    that the variables in the UGRID file are expected to have long names
    which are requirements over and above the UGRID specification.

    Parameters
    ----------

    filename : str
        Name of the UGRID mesh file to load.

    Returns
    -------

    : :class:`iris.cube.Cube`
        Cube has 1D coordinates for each of latitude and longitude, with the
        nodes from the original mesh translated into 4-element bounds.  The
        mesh topology, the face to node, and the face to face mapping are
        attached as attributes.  The data payload will be either an array of
        ones of the correct shape, or the data from the source file if a data
        constraint is provided.

    """
    if not isinstance(filename, str):
        if len(filename) != 1:
            raise ValueError(
                "Multiple mesh files provided to load_mesh: {}".format(filename)
            )
        filename = filename[0]

    # iris.load rather than ants.load to dodge derive_circular_status.
    # TODO:  revert back to ants.load and amend derive_circular_status to
    # ignore UGrid data (difficulty here is that we don't have a good test for
    # whether a cube is ugrid at this stage of the load).  Note: this will
    # also require updating patch in test_load_mesh.py
    cubes = iris.load(filename)
    ugrid = _UGridCubes(cubes)
    return ugrid._construct_mesh()


class _UGridCubes(object):
    """
    Extracts the information required for loading UGrid data from the source cubes.

    Parameters
    ----------
    cubes : :class:`iris.cube.CubeList`
        The cubes from which to extract the information.  Note that ANTS uses
        v2.3 of iris which uses ``extract`` with the ``strict`` argument
        rather than ``extract_cube``.
    data_constraint : :class:`iris.Constraint` (optional)
        Iris constraint to identify the data payload.  If omitted, cube will
        only contain the mesh and will not have a data payload.

    Attributes
    ----------
    cubes : :class:`iris.cube.CubeList`
        All the cubes from which the UGrid information is extracted.  Note
        that ANTS uses v2.3 of iris which uses ``extract`` with the ``strict``
        argument rather than ``extract_cube``.
    data_constraint : :class:`iris.Constraint` (optional)
        Iris constraint to identify the data payload.
    ugrid_cubes : :class:`iris.cube.CubeList`
        Those cubes identified as providing additional mesh information, as
        per UGrid spec.  This means the topology, face to node mapping and the
        face to face mapping.  Note that ANTS uses v2.3 of iris which uses
        ``extract`` with the ``strict`` argument rather than ``extract_cube``.

    """

    def __init__(self, cubes, data_constraint=None):
        self.cubes = cubes
        self.data_constraint = data_constraint
        self.ugrid_cubes = [
            cube
            for cube in cubes
            if cube.attributes.get("cf_role", "").lower() in UGRID_ROLES
        ]

    def _get_topologies(self):
        result = [
            cube
            for cube in self.cubes
            if cube.attributes.get("cf_role") == "mesh_topology"
        ]
        # TODO: Fix in #1185
        if len(result) != 1:
            raise RuntimeError(f"Expected exactly 1 topology, found {len(result)}.")
        return result

    @staticmethod
    def _get_face_coordinate_array(source, name):
        # Standard name ambiguous between face and nodes (e.g. 'latitude' is used
        # for both node latitude and face latitude).  Hence, reliance on
        # long_name() below
        try:
            # Implies source is a cubelist with individual lat/lon cubes
            result = source.extract(
                iris.Constraint(cube_func=lambda x: x.long_name == name), strict=True
            ).data
        except iris.exceptions.ConstraintMismatchError:
            # Implies source is a cubelist with regular cubes with lat/lon
            # coordinates
            coords = [
                c.coord(long_name=name) for c in source if c.coords(long_name=name)
            ]
            reference = coords[0]
            for coord in coords:
                if ants.utils.coord.relaxed_equality(reference, coord) is False:
                    raise ValueError(
                        "Trying to load a mesh with ambiguous coordinate {}".format(
                            name
                        )
                    )
            result = reference.points
        return result

    def _face_to_node_mapping_as_tuple(self, cube):
        """
        Derives and returns a namedtuple from the provided face to node mapping cube.

        Parameters
        ----------
        cube : :class:`iris.cube.Cube`
            Cube from which the face to node mapping is extracted.

        Returns
        -------
        : `collections.namedtuple`
            Contains two fields: `data` containing the face to node mapping
        data, and `metadata` containing metadata as a `dict`.  The metadata
        includes the cube attributes, var_name and long_name.

        """
        return _DataTuple(cube.data, self._make_dict(cube))

    @staticmethod
    def _make_dict(cube):
        # Need to convert iris LimitedAttributeDict to proper dict so we can set
        # long name (which, in turn, is needed on save).
        result = dict()
        result.update(cube.attributes)
        result["var_name"] = cube.var_name
        result["long_name"] = cube.long_name
        return result

    @staticmethod
    def get_nodes(node_longitudes, node_latitudes):
        """
        Derives and returns a namedtuple from the node locations

        Parameters
        ----------
        node_longitudes : :class:`np.ndarray`
            The longitudes of the nodes.
        node_latitudes : :class:`np.ndarray`
            The latitudes of the nodes.

        Returns
        -------
        : `collections.namedtuple`
            Contains two fields: `latitudes` and `longitudes` containing the
            respective locations for each node.

        """
        result = _NodeTuple(
            node_latitudes,
            node_longitudes,
        )
        return result

    @staticmethod
    def _construct_coordinate(nodes, face_centres, mapping, standard_name=None):
        """
        Construct a single iris coordinate from several cubes.

        Takes coordinate points from the face_centres data payload, and uses the
        nodes and mapping to generate the bounds.  Needs to be called separately
        per coordinate - i.e. once for latitude, and once for longitude.

        Parameters
        ----------

        nodes : :class:`iris.cube.Cube`
            The locations of the nodes from which to form the
            coordinate bounds (e.g. node_x for longitude coordinate).

        face_centres : :class:`iris.cube.Cube` or :class:`iris.coord.Coord`
            The locations of the face centres from which to form the
            coordinate points (e.g. face_x for longitude coordinate).  Note
            that ANTS uses v2.3 of iris which does not have the
            `nearest_neighbour_index` coordinate method.

        mapping : :class:`iris.cube.Cube`
            The mapping of the 4 nodes per face.

        standard_name : str (optional)
            CF Standard name for the coordinate.

        """
        node_positions = nodes.data
        bounds = node_positions[mapping.data - 1]

        units = nodes.units
        coord_system = iris.coord_systems.GeogCS(6371229.0)
        try:
            result = iris.coords.AuxCoord(
                face_centres.data, bounds=bounds, units=units, coord_system=coord_system
            )
        except AttributeError:
            """Face centres is a coordinate rather than a cube"""
            result = iris.coords.AuxCoord(
                face_centres.points,
                bounds=bounds,
                units=units,
                coord_system=coord_system,
            )
        if standard_name:
            result.standard_name = standard_name

        return result

    def _get_names(self, topology):
        names = {
            "latitude": {"standard": "latitude"},
            "longitude": {"standard": "longitude"},
        }
        node_var_names = topology["node_coordinates"].split()
        node_variable = [cube for cube in self.cubes if cube.var_name in node_var_names]
        node_lat_name = node_variable[0].long_name
        node_lon_name = node_variable[1].long_name
        if "lat" not in node_lat_name:
            # node_lat_name and node_lon_name could be in either order - if
            # they're wrong, let's correct them:
            node_lat_name, node_lon_name = node_lon_name, node_lat_name

        names["latitude"]["nodes"] = node_lat_name
        names["longitude"]["nodes"] = node_lon_name

        face_var_names = topology["face_coordinates"].split()
        face_variable = [cube for cube in self.cubes if cube.var_name in face_var_names]
        if len(face_variable) > 0:
            # Implies cubes for face coordinates (alternative is coordinates)
            face_lat_name = face_variable[0].long_name
            face_lon_name = face_variable[1].long_name
            if "lat" not in face_lat_name:
                # face_lat_name and face_lon_name could be in either order - if
                # they're wrong, let's correct them:
                face_lat_name, face_lon_name = face_lon_name, face_lat_name

            names["latitude"]["faces"] = face_lat_name
            names["longitude"]["faces"] = face_lon_name

        return names

    def _construct_mesh(self):
        topology = self._get_topologies()[0]
        names = self._get_names(topology.attributes)
        node_latitudes = self.cubes.extract(
            iris.Constraint(
                cube_func=lambda x: x.long_name == names["latitude"]["nodes"]
            ),
            strict=True,
        )
        node_longitudes = self.cubes.extract(
            iris.Constraint(
                cube_func=lambda x: x.long_name == names["longitude"]["nodes"]
            ),
            strict=True,
        )
        if self.data_constraint:
            # Face lats/lons are attached to data as coordinates rather than being
            # exposed as cubes in protocubes:
            result = self.cubes.extract(self.data_constraint, strict=True)
            face_latitudes = result.coord(axis="y")
            face_longitudes = result.coord(axis="x")
            # Why remove the existing lat/lon coordinates?  Ensures consistency
            # with data loaded from mesh, by generating coordinates from the same
            # face/node/mapping information (also ensures error if either node or
            # mapping information is missing):
            result.remove_coord(face_latitudes)
            result.remove_coord(face_longitudes)
            names["latitude"]["faces"] = face_latitudes.long_name
            names["longitude"]["faces"] = face_longitudes.long_name
            face_latitudes = face_latitudes.points
            face_longitudes = face_longitudes.points
        else:
            if "faces" not in names["latitude"]:
                # Implies face coordinates have loaded as coordinates
                # (i.e. there is at least one cube with data), and we're not
                # constraining the load to a particular data cube, so need to
                # derive the face latitudes from the entire cube list.
                latitude_names = {
                    cube.coord(axis="y").long_name
                    for cube in self.cubes
                    if len(cube.coords(axis="y")) == 1
                }
                if len(latitude_names) != 1:
                    raise RuntimeError(
                        f"Inconsistent latitude names {latitude_names} found during "
                        "load.  Try loading with a data constraint to pick out the "
                        "data cube of interest."
                    )
                names["latitude"]["faces"] = latitude_names.pop()
            if "faces" not in names["longitude"]:
                # Implies face coordinates have loaded as coordinates
                # (i.e. there is at least one cube with data), and we're not
                # constraining the load to a particular data cube, so need to
                # derive the face longitudes from the entire cube list.
                longitude_names = {
                    cube.coord(axis="x").long_name
                    for cube in self.cubes
                    if len(cube.coords(axis="x")) == 1
                }
                if len(longitude_names) != 1:
                    raise RuntimeError(
                        f"Inconsistent longitude names {longitude_names} found during "
                        "load.  Try loading with a data constraint to pick out the "
                        "data cube of interest."
                    )
                names["longitude"]["faces"] = longitude_names.pop()
            face_latitudes = self._get_face_coordinate_array(
                self.cubes, names["latitude"]["faces"]
            )
            face_longitudes = self._get_face_coordinate_array(
                self.cubes, names["longitude"]["faces"]
            )
            _result_shape = face_latitudes.shape
            result = iris.cube.Cube(data=da.ones(_result_shape))

        face_to_node_mapping = self.cubes.extract(
            iris.Constraint(
                cube_func=lambda x: "cf_role" in x.attributes
                and x.attributes["cf_role"] == "face_node_connectivity"
            ),
            strict=True,
        )
        latitude = self._construct_coordinate(
            node_latitudes,
            face_latitudes,
            face_to_node_mapping,
            names["latitude"]["standard"],
        )
        longitude = self._construct_coordinate(
            node_longitudes,
            face_longitudes,
            face_to_node_mapping,
            names["longitude"]["standard"],
        )

        latitude.long_name = names["latitude"]["faces"]
        longitude.long_name = names["longitude"]["faces"]
        latitude.attributes["bounds_long_name"] = names["latitude"]["nodes"]
        longitude.attributes["bounds_long_name"] = names["longitude"]["nodes"]

        result.add_aux_coord(latitude, len(result.shape) - 1)
        result.add_aux_coord(longitude, len(result.shape) - 1)

        # Immutability isn't enforced - but none of our processing operations
        # should touch this metadata.  It's needed on save, to avoid needing to
        # recompute from the coordinates etc
        immutable_metadata = {}
        immutable_metadata[topology.attributes["cf_role"]] = self._make_dict(topology)
        immutable_metadata[
            face_to_node_mapping.attributes["cf_role"]
        ] = self._face_to_node_mapping_as_tuple(face_to_node_mapping)
        immutable_metadata["nodes"] = self.get_nodes(
            node_latitudes=node_latitudes.data, node_longitudes=node_longitudes.data
        )
        result.attributes.update(immutable_metadata)

        # Should we include a CF version here as well?
        result.attributes["Conventions"] = "UGRID-1.0"

        return result


class _DataTuple(namedtuple("Mapping", "data metadata")):
    def __eq__(self, other):
        for field in self._fields:
            self_field = getattr(self, field)
            try:
                other_field = getattr(other, field)
            except AttributeError:
                # Missing attribute in other - rather than exception, makes
                # sense to consider other to be not equal to self:
                return False
            if isinstance(self_field, np.ndarray):
                equal = np.allclose(self_field, other_field)
            else:
                equal = self_field == other_field
            if equal is False:
                # Early return to shortcircuit once we hit any unequal field.
                return equal
        return True

    def __ne__(self, other):
        return not (self == other)


class _NodeTuple(namedtuple("Nodes", "latitudes longitudes")):
    def __eq__(self, other):
        for field in self._fields:
            self_field = getattr(self, field)
            try:
                other_field = getattr(other, field)
            except AttributeError:
                # Missing attribute in other - rather than exception, makes
                # sense to consider other to be not equal to self:
                return False
            if isinstance(self_field, np.ndarray):
                equal = np.allclose(self_field, other_field)
            else:
                equal = self_field == other_field
            if equal is False:
                # Early return to shortcircuit once we hit any unequal field.
                return equal
        return True

    def __ne__(self, other):
        return not (self == other)


def _get_lat_lon(cube, face_lats_variable, face_lons_variable):
    """
    Gets coordinates in format suitable for saving.

    In particular, this means that the bounds are removed (UGrid nodes are not
    created from the bounds yet.  Instead, they are saved from the information
    in the cube nodes attribute).

    Parameters
    ----------
    cube: :class:`iris.cube.Cube`
        Cube from which to extract the coordinates.
    face_lats_variable: str
        Name of the latitude variable.
    face_lons_variable: str
        Name of the longitude variable.

    Returns
    -------
    : tuple of (:class:`iris.coord.Coord`, :class:`iris.coord.Coord`)
        Latitude and longitude coordinates in a tuple, with the bounds removed.
    """
    latitude = cube.coord(face_lats_variable)
    longitude = cube.coord(face_lons_variable)

    return latitude, longitude


def _valid_cube_names(cubes):
    """
    Makes cube names valid if possible, or rejects if not.

    Purpose of this is to avoid duplicate cube names (iris creates unique
    names on save, but these are not fed back to cube, so it's hard to
    identify them - but not impossible, if we have a compelling need to do
    so).  So, for current UGrid purposes, we're defining valid as unique
    names.

    Operates on cubes in place.

    Parameters
    ----------
    cubes: :class:`iris.cube.CubeList`
        List of cubes with names to be validated.  Note that ANTS uses v2.3 of
        iris which uses ``extract`` with the ``strict`` argument rather than
        ``extract_cube``.

    """
    # If our name translation above lets duplicate names through, let's throw
    # an error (will need to relax this constraint later, most likely).
    names = {cube.name() for cube in cubes}
    if len(names) != len(cubes):
        raise IOError(
            "Cannot save cubes with same name.  Saw cubes with "
            "names: {}".format(names)
        )


def _save_nodes(dataset, variable_name, faces, nodes, units):
    """
    Saves the node information.

    Parameters
    ----------
    dataset : :class:`netCDF4.dataset`
        Dataset to save the node information to.
    variable_name : str
        name of :class:`netCDF4.variable` to save the node information to.
    faces : :class:`iris.coords.Coord`
        Coordinate from which to take the metadata for the node -
        typically the corresponding face coordinate (hence parameter
        name).
    nodes : :class:`numpy.ndarray`
        Array with the node values for the coordinate being saved.
    units : str
        Name for the units - can't be inferred from faces parameter.

    Returns
    -------
    : None
    """

    node_variable = dataset.createVariable(
        variable_name, np.float64, dimensions=("num_node",)
    )
    node_variable.standard_name = faces.standard_name
    node_variable.long_name = faces.attributes["bounds_long_name"]
    node_variable[:] = nodes
    node_variable.units = units


def _save_connectivity(dataset, cube, face_node_variable, face_dim):
    """
    Saves the face to node mapping information.

    Parameters
    ----------
    dataset : :class:`netCDF4.dataset`
        Dataset to save the face to node mapping to.
    cube : :class:`iris.cube.Cube`
        UGrid format cube to save.
    face_node_variable : str
        Name of the netCDF variable to save the face node mapping information
        to.
    face_dim : int
        Index of the dimension to save the face_node_variable against.

    Returns
    -------
    : None

    """
    # np.int32 is derived from CF allowed types:
    # http://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/cf-conventions.html#_data_types
    # which in turn leads to
    # https://www.unidata.ucar.edu/software/netcdf/docs/netcdf_8h.html#a306f8e52ab0aae4b5b902ada169b7b3c
    # Notice that CF only includes 16 and 32 bit ints, and we'd exhaust 16bit
    # ints at pretty moderate resolutions (C100 mesh has 6x100x100 faces, so
    # is in the ball park for the 16bit int limit).
    face_node_mapping = dataset.createVariable(
        face_node_variable, np.int32, dimensions=(face_dim, "num_vertices")
    )
    metadata = cube.attributes["face_node_connectivity"].metadata
    face_node_mapping.setncatts(metadata)

    face_node_mapping[:] = cube.attributes["face_node_connectivity"].data


def _save_minimal_topology(dataset, name):
    """
    Saves just enough topology for LFRic load to work.

    Parameters
    ----------
    dataset : :class:`netCDF4.dataset`
        Dataset to save the topology to.
    name : str
        Name of the netCDF variable to save the topology to.

    Returns
    -------
    : :class:`netCDF4.Variable`
        The topology variable.  Can be amended until the dataset is closed.

    """
    # np.int32 is derived from CF allowed types:
    # http://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/cf-conventions.html#_data_types
    # which in turn leads to
    # https://www.unidata.ucar.edu/software/netcdf/docs/netcdf_8h.html#a306f8e52ab0aae4b5b902ada169b7b3c
    # Notice that CF only includes 16 and 32 bit ints, and we'd exhaust 16bit
    # ints at pretty moderate resolutions (C100 mesh has 6x100x100 faces, so
    # is in the ball park for the 16bit int limit).
    topology = dataset.createVariable(name, np.int32, dimensions=())
    topology.cf_role = "mesh_topology"

    return topology


def _save_UGRID_compliance(dataset, cubes, topology, coordinates, name):
    # Additional metadata needed for CF/UGRID compliance:
    ref_cube = cubes[0]
    dataset.Conventions = ref_cube.attributes["Conventions"]
    topology.topology_dimension = 2
    topology.node_coordinates = "{} {}".format(
        coordinates["nodes"]["latitudes"], coordinates["nodes"]["longitudes"]
    )
    topology.face_coordinates = "{} {}".format(
        coordinates["faces"]["latitudes"], coordinates["faces"]["longitudes"]
    )
    topology.face_node_connectivity = coordinates["face_node_mapping"]["name"]
    topology.face_dimension = coordinates["faces"]["dim"]
    topology.long_name = ref_cube.attributes["mesh_topology"]["long_name"]


def check_validity(cubes, face_lats_variable, face_lons_variable):
    """
    Checks that the conditions to be able to save the UGrid cube are met.

    At present, these checks are:

    1. No saver argument was provided.  This raises a warning if a saver
    argument was provided.

    2. All cubes have the same UGrid attributes: topology, face to node
    mapping, nodes and conventions.  This raises an error if the attributes
    differ.

    3. All cubes have the same latitude and longitude coordinates (to within
    floating point tolerances).  This raises an error if the coordinates
    differ.

    Parameters
    ----------
    cubes: :class:`iris.cube.CubeList`
        Cubes to validate.  Note that ANTS uses v2.3 of iris which uses
        ``extract`` with the ``strict`` argument rather than ``extract_cube``.
    face_lats_variable: str
        Name of the latitude variable.
    face_lons_variable: str
        Name of the longitude variable.

    Returns
    -------

    : None
        Operates in place.

    """
    _valid_cube_names(cubes)
    ref_cube = cubes[0]
    # 1st docstring criterion: warning if a saver is given
    if ants.config.CONFIG["saver"] is not None:
        msg = (
            "Invalid fileformat {} specified.  Saving UGrid cubes {} as "
            "UGrid format.".format(
                ants.config.CONFIG["saver"], " ".join([cube.name() for cube in cubes])
            )
        )
        warnings.warn(msg)
    # 2nd docstring criterion: error if UGrid attributes differ
    expected_attributes = (
        "mesh_topology",
        "face_node_connectivity",
        "Conventions",
        "nodes",
    )
    for expected_attribute in expected_attributes:
        reference = ref_cube.attributes[expected_attribute]
        for cube in cubes:
            if cube.attributes[expected_attribute] != reference:
                raise IOError(
                    "Cannot save cubes with different {} attribute.".format(
                        expected_attribute
                    )
                )
    # 3rd docstring criterion: error if lat or lon coords differ
    if ants.utils.cube.is_equal_hgrid(cubes) is False:
        msg = (
            "Cannot save multiple UGrid cubes with different latitude "
            "or longitude coordinates."
        )
        raise ValueError(msg)


def _apply_XIOS_workaround(cubes):
    """
    XIOS requires a long name for each phenomenon.

    Parameters
    ----------

    cubes : :class:`iris.cube.CubeList`
        UGrid cubes to be fixed in-place for XIOS compliance.  Note that ANTS
        uses v2.3 of iris which uses ``extract`` with the ``strict`` argument
        rather than ``extract_cube``.

    """
    for cube in cubes:
        # XIOS only cares that there is a long name, not what the long name is:
        if cube.long_name is None:
            cube.long_name = cube.standard_name
        # The value of online operation is important though.  There are other
        # valid values, so let the user manually set one on the cube before
        # saving.  If the user has manually set a value, let's assume they
        # know why they've chosen the value and not warn them that they're not
        # getting the default:
        if "online_operation" not in cube.attributes:
            cube.attributes["online_operation"] = "once"
        # Value of coordinate var_names may or may not be significant.
        # Because of this uncertainty, it may be that the default value is the
        # only valid value.  So let's warn the user if the default isn't used
        # (if we later find out what the limits of XIOS actually are, we can
        # revise this to clobber the var_name or ignore it as appropriate):
        for axis in ("x", "y"):
            coord = cube.coord(axis=axis)
            default_var_name = (
                f"{cube.attributes['mesh_topology']['var_name']}_face_{axis}"
            )
            if coord.var_name is None:
                coord.var_name = default_var_name
            elif coord.var_name != default_var_name:
                msg = warnings.UserWarning(
                    f"Coordinate {coord.name()} has an existing var_name of "
                    "{coord.var_name}.  This has been left intact rather than being "
                    "replaced with a known good var_name.  To get the default "
                    "coordinate var_name, remove the var_name prior to saving."
                )
                warnings.warn(msg)


def iris_save(cubes, filename, face_lats_variable, face_lons_variable, mesh_name):
    """
    Contains the parts of the save that are handled by iris.

    cube : :class:`iris.cube.Cube` or :class:`iris.cube.CubeList`
        UGrid cube or cubes with the data to be saved.  Note that ANTS uses
        v2.3 of iris which uses ``extract`` with the ``strict`` argument
        rather than ``extract_cube``.

    filename : str
        Filename for the UGRID mesh file.

    """
    # Need to alter cube to remove the parts iris can't currently save, so to
    # prevent altering original, we'll take a copy of the cube.  Current LFRic
    # constraint is that data type must be float64:
    # TODO: Support different datatypes
    cubes = iris.cube.CubeList([cube.copy(np.float64(cube.data)) for cube in cubes])

    for cube in cubes:
        # Drop attributes defined as UGRID_ROLES since these aren't handled by
        # iris as we need them to be.  Instead, these are handled in
        # _save_connectivity and _save_minimal_topology
        [cube.attributes.pop(role) for role in UGRID_ROLES if role in cube.attributes]

        # But we need to add the UGrid spec attributes for describing the data
        # variables:
        cube.attributes["location"] = "face"
        cube.attributes["mesh"] = mesh_name

        # Drop face_node_connectivity attribute - it's already saved as a
        # Variable in _save_connectivity:
        if "face_node_connectivity" in cube.attributes:
            cube.attributes.pop("face_node_connectivity")

        # Drop node attributes, handled separately:
        cube.attributes.pop("nodes")

        # Can safely remove the lat/lon bounds (nodes are saved from the
        # attributes on the cubes, so the bounds information just pollutes the
        # save):
        latitude, longitude = _get_lat_lon(cube, face_lats_variable, face_lons_variable)
        latitude.bounds = longitude.bounds = None
        # Drop bounds_long_name attribute from coordinates: it's only purpose
        # is in our additional save to link the nodes to the correct faces,
        # but it causes issues with XIOS if it persists to disc.
        [
            coord.attributes.pop("bounds_long_name")
            for coord in (latitude, longitude)
            if "bounds_long_name" in coord.attributes
        ]
        # Drop coordinate reference systems (not needed yet; want to keep
        # nodes and faces consistent, so simpler to remove crs for faces from
        # iris save than to add crs for nodes to additional save).
        latitude.coord_system = longitude.coord_system = None

    iris.save(
        cubes,
        filename,
        saver="nc",
        local_keys=("location", "mesh", "online_operation"),
        unlimited_dimensions=[],
        # Ticket #1829: Fill value should be largest finite negative value.
        fill_value=np.nan_to_num(np.NINF),
    )


def additional_save(cubes, filename, face_lats_variable, face_lons_variable, mesh_name):
    """
    Contains the parts of the save that are not handled by iris.

    Parameters
    ----------
    cube : :class:`iris.cube.Cube`
        UGrid cube with the data to be saved.
    filename : str
        Filename for the UGRID mesh file.
    face_lats_variable: str
        Name of the latitude variable.
    face_lons_variable: str
        Name of the longitude variable.

    Returns
    -------
    : None
    """
    ref_cube = cubes[0]
    latitude, longitude = _get_lat_lon(ref_cube, face_lats_variable, face_lons_variable)
    dataset = Dataset(filename, "a")

    topology = _save_minimal_topology(dataset, mesh_name)

    nodes = ref_cube.attributes["nodes"]
    dataset.createDimension("num_node", len(nodes.latitudes))
    dataset.createDimension("num_vertices", 4)

    # ugrid.node_lat does not contain a unit - instead, we'll hard-code CF
    # unit (doesn't seem to be an easy way to get 'degrees_north' out of
    # the face coordinate):
    node_lats_variable = f"{topology.name}_node_y"
    _save_nodes(dataset, node_lats_variable, latitude, nodes.latitudes, "degrees_north")

    node_lons_variable = f"{topology.name}_node_x"
    _save_nodes(
        dataset, node_lons_variable, longitude, nodes.longitudes, "degrees_east"
    )
    face_node_name = ref_cube.attributes["face_node_connectivity"].metadata.pop(
        "var_name"
    )
    face_dim = "dim{}".format(ref_cube.coord_dims("latitude")[0])
    _save_connectivity(dataset, ref_cube, face_node_name, face_dim)

    coordinate_info = {
        "nodes": {"latitudes": node_lats_variable, "longitudes": node_lons_variable},
        "faces": {
            "latitudes": latitude.var_name,
            "longitudes": longitude.var_name,
            "dim": face_dim,
        },
        "face_node_mapping": {"name": face_node_name},
    }

    _save_UGRID_compliance(dataset, cubes, topology, coordinate_info, mesh_name)

    dataset.close()


def save(cubes, filename):
    """
    Saves a cube or a cubelist as a ugrid mesh file.

    Data payloads are from the cubes directly, but the coordinates rely on the
    cubes having node coordinates as attributes (guaranteed if the cubes
    were created with :func:`ants.fileformats._ugrid.load_mesh`).

    It's the users responsibility to ensure their processing has not caused
    the cubes' metadata and coordinates to become out of sync to those on the
    node attributes, and to rectify any discrepancies.

    Parameters
    ----------

    cube : :class:`iris.cube.Cube` or :class:`iris.cube.CubeList`
        UGrid cube or cubes with the data to be saved.  Note that ANTS uses
        v2.3 of iris which uses ``extract`` with the ``strict`` argument
        rather than ``extract_cube``.

    filename : str
        Filename for the UGRID mesh file.

    """
    cubes = ants.utils.cube.as_cubelist(cubes)

    # Why create these variables?  Need them to ensure consistency between
    # netCDF variable name created by iris, and reference to the variable in
    # the topology
    face_lats_variable = "latitude"
    face_lons_variable = "longitude"

    # Dev note: The iris part of the save handles cube lists.  Most of the
    # extra part of the save, using netCDF4 python directly, does not, and
    # only uses the first cube.  However, all information in this part of the
    # save is identical for all cubes (verified in check_validity).  If
    # additional_save is amended, the validity check should also be updated to
    # ensure it's consistent with the additional_save behaviour.
    check_validity(cubes, face_lats_variable, face_lons_variable)
    mesh_name = cubes[0].attributes["mesh_topology"]["var_name"]

    _apply_XIOS_workaround(cubes)

    iris_save(cubes, filename, face_lats_variable, face_lons_variable, mesh_name)
    additional_save(cubes, filename, face_lats_variable, face_lons_variable, mesh_name)
