# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import unittest.mock as mock
from tempfile import NamedTemporaryFile

import ants.tests
import cf_units
import iris
import numpy as np
from ants.fileformats._ugrid import _DataTuple, _NodeTuple, load, load_mesh, save
from netCDF4 import Dataset

TOPOLOGY_NAME = ants.tests.stock.mesh_C4().attributes["mesh_topology"]["var_name"]


class _LoadMeshCommon(object):
    # Loading mesh from a mesh file w/o a data payload.
    def setUp(self):
        self.cube = load_mesh(ants.tests.get_data_path(("stock", "mesh_C4.nc")))


class _LoadDataCommon(object):
    # Loading mesh and data payloads from a file with a data payload.
    def setUp(self):
        # In this case, self.cube is a CubeList, not a Cube.  Need self.cube
        # name for compatibility with other mixin classes.
        self.cube = load(ants.tests.get_data_path(("stock", "data_C4.nc")))


class _LoadMeshFromDataCommon(object):
    # Loading mesh only from a file with a data payload.
    def setUp(self):
        self.cube = load_mesh(ants.tests.get_data_path(("stock", "data_C4.nc")))


class _CoordinateLoadCommon(object):
    def test_points_shape(self):
        expected = (96,)

        for cube in ants.utils.cube.as_cubelist(self.cube):
            coordinate = cube.coord(self.coord_name)
            with self.subTest():
                self.assertEqual(coordinate.points.shape, expected)

    def test_bounds_shape(self):
        expected = (96, 4)

        for cube in ants.utils.cube.as_cubelist(self.cube):
            coordinate = cube.coord(self.coord_name)
            with self.subTest():
                self.assertEqual(coordinate.bounds.shape, expected)

    def test_bounds_range(self):
        expected = self.bounds_limits

        for cube in ants.utils.cube.as_cubelist(self.cube):
            coordinate = cube.coord(self.coord_name)
            actual = (np.min(coordinate.bounds), np.max(coordinate.bounds))

            with self.subTest():
                self.assertArrayAlmostEqual(actual, expected)

    def test_points_range(self):
        expected = self.points_limits

        for cube in ants.utils.cube.as_cubelist(self.cube):
            coordinate = cube.coord(self.coord_name)
            actual = (np.min(coordinate.points), np.max(coordinate.points))

            with self.subTest():
                self.assertArrayAlmostEqual(actual, expected)

    def test_units(self):
        expected = self.units

        for cube in ants.utils.cube.as_cubelist(self.cube):
            coordinate = cube.coord(self.coord_name)

            with self.subTest():
                self.assertEqual(coordinate.units, expected)

    def test_coordinate_system(self):
        expected = self.coord_system

        for cube in ants.utils.cube.as_cubelist(self.cube):
            coordinate = cube.coord(self.coord_name)

            with self.subTest():
                self.assertEqual(coordinate.coord_system, expected)

    def test_nodes_name(self):
        expected = "{} of 2D mesh nodes.".format(self.coord_name)

        for cube in ants.utils.cube.as_cubelist(self.cube):
            coordinate = cube.coord(self.coord_name)
        actual = coordinate.attributes["bounds_long_name"]

        with self.subTest():
            self.assertEqual(actual, expected)

    def test_long_name(self):
        expected = "{} of 2D face centres".format(self.coord_name)

        for cube in ants.utils.cube.as_cubelist(self.cube):
            coordinate = cube.coord(self.coord_name)
        actual = coordinate.long_name

        with self.subTest():
            self.assertEqual(actual, expected)


class TestLoadMeshLatitude(_CoordinateLoadCommon, _LoadMeshCommon, ants.tests.TestCase):

    coord_name = "latitude"
    points_limits = (-74.247943, 74.247943)
    bounds_limits = (-90.0, 90.0)
    units = cf_units.Unit("degrees_north")
    coord_system = iris.coord_systems.GeogCS(6371229.0)


class TestLoadMeshLongitude(
    _CoordinateLoadCommon, _LoadMeshCommon, ants.tests.TestCase
):

    coord_name = "longitude"
    # 16 cells around equatorial panels (4x4 for C4_mesh); so max
    # longitude is 360. - (360./16.) = 337.5.  In other words, nodes on
    # meridian defined as 0 degrees (not 360), so cells that abut meridian
    # from the West have two nodes at ~337.5 degrees and two at 0 degrees.
    bounds_limits = (0.0, 337.5)
    points_limits = (11.279747, 348.720253)
    units = cf_units.Unit("degrees_east")
    coord_system = iris.coord_systems.GeogCS(6371229.0)


class TestLoadDataLatitude(_CoordinateLoadCommon, _LoadDataCommon, ants.tests.TestCase):

    coord_name = "latitude"
    points_limits = (-74.247943, 74.247943)
    bounds_limits = (-90.0, 90.0)
    units = cf_units.Unit("degrees_north")
    coord_system = iris.coord_systems.GeogCS(6371229.0)


class TestLoadDataLongitude(
    _CoordinateLoadCommon, _LoadDataCommon, ants.tests.TestCase
):

    coord_name = "longitude"
    # 16 cells around equatorial panels (4x4 for C4_mesh); so max
    # longitude is 360. - (360./16.) = 337.5.  In other words, nodes on
    # meridian defined as 0 degrees (not 360), so cells that abut meridian
    # from the West have two nodes at ~337.5 degrees and two at 0 degrees.
    bounds_limits = (0.0, 337.5)
    points_limits = (11.279747, 348.720253)
    units = cf_units.Unit("degrees_east")
    coord_system = iris.coord_systems.GeogCS(6371229.0)


class TestLoadMeshFromDataLatitude(
    _CoordinateLoadCommon, _LoadMeshFromDataCommon, ants.tests.TestCase
):

    coord_name = "latitude"
    points_limits = (-74.247943, 74.247943)
    bounds_limits = (-90.0, 90.0)
    units = cf_units.Unit("degrees_north")
    coord_system = iris.coord_systems.GeogCS(6371229.0)


class TestLoadMeshFromDataLongitude(
    _CoordinateLoadCommon, _LoadMeshFromDataCommon, ants.tests.TestCase
):

    coord_name = "longitude"
    # 16 cells around equatorial panels (4x4 for C4_mesh); so max
    # longitude is 360. - (360./16.) = 337.5.  In other words, nodes on
    # meridian defined as 0 degrees (not 360), so cells that abut meridian
    # from the West have two nodes at ~337.5 degrees and two at 0 degrees.
    bounds_limits = (0.0, 337.5)
    points_limits = (11.279747, 348.720253)
    units = cf_units.Unit("degrees_east")
    coord_system = iris.coord_systems.GeogCS(6371229.0)


class TestLoadMeshMetaData(_LoadMeshCommon, ants.tests.TestCase):
    def test_face_node_mapping_loaded_as_attribute(self):
        self.assertIn("face_node_connectivity", self.cube.attributes)

    def test_face_node_mapping_type(self):
        actual = self.cube.attributes["face_node_connectivity"]
        self.assertIsInstance(actual, ants.fileformats._ugrid._DataTuple)

    def test_face_node_mapping_data(self):
        dataset = Dataset(ants.tests.get_data_path(("stock", "mesh_C4.nc")))
        # If stock mesh is updated, this hard-coded name may need updating too:
        expected = dataset[f"{TOPOLOGY_NAME}_face_nodes"][:]
        actual = self.cube.attributes["face_node_connectivity"].data

        self.assertArrayEqual(expected, actual)

    def test_face_node_mapping_data_type(self):
        actual = self.cube.attributes["face_node_connectivity"].data
        self.assertIsInstance(actual, np.ndarray)

    def test_face_node_mapping_metadata_type(self):
        actual = self.cube.attributes["face_node_connectivity"].metadata
        self.assertIsInstance(actual, dict)

    def test_node_longitudes(self):
        dataset = Dataset(ants.tests.get_data_path(("stock", "mesh_C4.nc")))
        # If stock mesh is updated, this hard-coded name may need updating too:
        expected = dataset[f"{TOPOLOGY_NAME}_node_x"][:]

        actual = self.cube.attributes["nodes"].longitudes

        self.assertArrayEqual(actual, expected)

    def test_node_latitudes(self):
        dataset = Dataset(ants.tests.get_data_path(("stock", "mesh_C4.nc")))
        # If stock mesh is updated, this hard-coded name may need updating too:
        expected = dataset[f"{TOPOLOGY_NAME}_node_y"][:]

        actual = self.cube.attributes["nodes"].latitudes

        self.assertArrayEqual(actual, expected)

    def test_conventions(self):
        expected = "UGRID-1.0"

        actual = self.cube.attributes["Conventions"]

        self.assertEqual(actual, expected)

    def test_variable_name(self):
        expected = "example_C4"

        actual = self.cube.attributes["mesh_topology"]["var_name"]

        self.assertEqual(actual, expected)


class TestLoadDataMetaData(ants.tests.TestCase):
    def setUp(self):
        self.cubes = load(ants.tests.get_data_path(("stock", "data_C4.nc")))

    def _assert_attribute(self, name):
        # Slight nuance here - this checks multiple things with a single
        # assertion for each cube in self.cubes:
        # 1. the name attribute exists,
        # 2. the name attribute is a dict (or at least has a
        # __getitem__).
        # 3. the name of the attribute on the cube is also the cf_role item in
        # the dict.  Given cf_roles are constrained by spec, makes sense to
        # ensure the cube attribute name also follows that spec (until there's
        # a reason to relax this constraint).
        for cube in self.cubes:
            with self.subTest():
                self.assertEqual(name, cube.attributes[name]["cf_role"])

    def test_returns_cubelist(self):
        expected = iris.cube.CubeList

        actual = self.cubes

        self.assertIsInstance(actual, expected)

    def test_mesh_topology_loaded_as_attribute(self):
        self._assert_attribute("mesh_topology")

    def test_face_face_mapping_not_loaded_as_attribute(self):
        # Unlike a mesh, there should be no face to face connectivity for a
        # UGrid data file (we haven't yet implemented a save for face to face
        # connectivity - this test reminds us to consider implementing it for
        # load after we implement for save).
        expected = "face_face_connectivity"

        for cube in self.cubes:
            actual = cube.attributes
            with self.subTest():
                self.assertNotIn(expected, actual)

    def test_face_node_mapping_loaded_as_attribute(self):
        for cube in self.cubes:
            actual = cube.attributes
            with self.subTest():
                self.assertIn("face_node_connectivity", actual)

    def test_face_node_mapping_type(self):
        for cube in self.cubes:
            actual = cube.attributes
            with self.subTest():
                actual = cube.attributes["face_node_connectivity"]
                self.assertIsInstance(actual, ants.fileformats._ugrid._DataTuple)

    def test_face_node_mapping_data(self):
        dataset = Dataset(ants.tests.get_data_path(("stock", "mesh_C4.nc")))
        # If stock mesh is updated, this hard-coded name may need updating too:
        expected = dataset[f"{TOPOLOGY_NAME}_face_nodes"][:]

        for cube in self.cubes:
            actual = cube.attributes["face_node_connectivity"].data
            with self.subTest():
                self.assertArrayEqual(expected, actual)

    def test_face_node_mapping_data_type(self):
        for cube in self.cubes:
            actual = cube.attributes["face_node_connectivity"].data
            with self.subTest():
                self.assertIsInstance(actual, np.ndarray)

    def test_face_node_mapping_metadata_type(self):
        for cube in self.cubes:
            actual = cube.attributes["face_node_connectivity"].metadata
            with self.subTest():
                self.assertIsInstance(actual, dict)

    def test_node_longitudes(self):
        dataset = Dataset(ants.tests.get_data_path(("stock", "data_C4.nc")))
        # If stock mesh is updated, this hard-coded name may need updating too:
        expected = dataset["example_C4_node_x"][:]

        for cube in self.cubes:
            actual = cube.attributes["nodes"].longitudes
            with self.subTest():
                self.assertArrayEqual(actual, expected)

    def test_node_latitudes(self):
        dataset = Dataset(ants.tests.get_data_path(("stock", "data_C4.nc")))
        # If stock mesh is updated, this hard-coded name may need updating too:
        expected = dataset["example_C4_node_y"][:]

        for cube in self.cubes:
            actual = cube.attributes["nodes"].latitudes
            with self.subTest():
                self.assertArrayEqual(actual, expected)

    def test_conventions(self):
        expected = "UGRID-1.0"

        for cube in self.cubes:
            actual = cube.attributes["Conventions"]
            with self.subTest():
                self.assertEqual(actual, expected)


class TestLoadData(ants.tests.TestCase):
    def test_mesh_only_file_without_constraint(self):
        expected = np.ones(96)

        actual = ants.tests.stock.mesh_C4(load_data=False).data

        self.assertArrayEqual(actual, expected)

    def test_data_file_loads_two_cubes(self):
        expected = 2

        actual = len(ants.tests.stock.mesh_C4(load_data=True, constraint=None))

        self.assertEqual(actual, expected)

    # Following two tests load in all data variables from our UGrid sample
    # data, and check that the relevant data array corresponds to the cube of
    # the correct name.
    def test_data_in_sample_data_without_constraint(self):
        expected_data = np.arange(1, 97)

        actual_cubes = ants.tests.stock.mesh_C4(load_data=True, constraint=None)
        actual_data = actual_cubes.extract_strict("sample_data").data

        self.assertArrayEqual(actual_data, expected_data)

    def test_data_in_additional_sample_data_without_constraint(self):
        expected_data = np.arange(1, 97) * -1

        actual_cubes = ants.tests.stock.mesh_C4(load_data=True, constraint=None)
        actual_data = actual_cubes.extract_strict("additional_sample_data").data

        self.assertArrayEqual(actual_data, expected_data)

    def test_return_type_without_constraint(self):
        expected = iris.cube.CubeList

        actual = ants.tests.stock.mesh_C4(load_data=True)

        self.assertIsInstance(actual, expected)

    def test_data_in_sample_data_with_constraint(self):
        expected_data = np.arange(1, 97)

        actual_cubes = ants.tests.stock.mesh_C4(
            load_data=True, constraint="sample_data"
        )
        self.assertEqual(len(actual_cubes), 1)
        actual_data = actual_cubes[0].data

        self.assertArrayEqual(actual_data, expected_data)

    def test_data_in_additional_sample_data_with_constraint(self):
        expected_data = np.arange(1, 97) * -1

        actual_cubes = ants.tests.stock.mesh_C4(
            load_data=True, constraint="additional_sample_data"
        )
        self.assertEqual(len(actual_cubes), 1)
        actual_data = actual_cubes[0].data

        self.assertArrayEqual(actual_data, expected_data)

    def test_multiple_constraints(self):
        expected_names = sorted(["sample_data", "additional_sample_data"])

        actual_cubes = ants.tests.stock.mesh_C4(
            load_data=True, constraint=["additional_sample_data", "sample_data"]
        )
        actual_names = sorted([c.name() for c in actual_cubes])

        self.assertEqual(expected_names, actual_names)

    def test_return_type_with_constraint(self):
        expected = iris.cube.CubeList

        actual = ants.tests.stock.mesh_C4(load_data=True, constraint="sample_data")

        self.assertIsInstance(actual, expected)

    def test_cube_can_be_copied(self):
        # Implication here is that there all attributes on cube can be
        # pickled.
        expected = iris.cube.Cube

        actual = ants.tests.stock.mesh_C4().copy()
        # Assertion here is really "assert doesn't raise", but since that
        # doesn't exist...
        self.assertIsInstance(actual, expected)


class _SaveCommon(object):
    def setUp(self):
        self.cube = ants.tests.stock.mesh_C4()
        # Need to store a reference to the payload so we can easily identify
        # it in later tests:
        self.data_name = "foo"

        self.cube.long_name = self.data_name
        # Doesn't matter what the dtype is, as long as it's not float64 (we
        # need to test later that data is converted to float64)
        self.cube.data = np.arange(1, 97, dtype=np.int32)
        self._temp_file = NamedTemporaryFile("w+b")
        save(self.cube, self._temp_file.name)
        self.result_dataset = Dataset(self._temp_file.name)

    def tearDown(self):
        self.result_dataset.close()
        self._temp_file.close()


class TestSaveData(_SaveCommon, ants.tests.TestCase):
    def test_data_values(self):
        expected = np.arange(1, 97)

        actual = self.result_dataset.variables[self.data_name][:]

        self.assertArrayEqual(expected, actual)

    def test_data_type(self):
        expected = np.float64

        actual = self.result_dataset.variables[self.data_name].dtype

        self.assertEqual(expected, actual)

    def test_has_long_name(self):
        expected = "long_name"

        actual = self.result_dataset.variables[self.data_name].ncattrs()

        self.assertIn(expected, actual)

    def test_has_correct_online_operation_attribute(self):
        # XIOS requirement, not part of CF/UGrid specs.
        expected = "once"
        actual = self.result_dataset.variables[self.data_name].online_operation

        self.assertEqual(expected, actual)


class _CoordinateSaveCommon(_SaveCommon):
    def _get_coordinate_variables(self, axis):
        # Only returns face and node variables.  Specifically, means the face
        # to node mapping is not present, so ensures there's no variable with
        # both 'face' and 'node' in the long name in subsequent tests.
        variables = self.result_dataset.variables
        return [
            variables[v]
            for v in variables
            if (
                "standard_name" in variables[v].ncattrs()
                and variables[v].standard_name.lower()
            )
            == axis
        ]

    def _get_variable(self, loci, axis):
        result = [
            v
            for v in self._get_coordinate_variables(axis)
            if loci in v.long_name.lower()
        ]
        if len(result) > 1:
            raise RuntimeError(f"Multiple variables found for {loci} on {axis} axis.")
        elif len(result) == 0:
            raise RuntimeError(f"No variables found for {loci} on {axis} axis.")
        return result[0]

    def test_no_bounds(self):
        # Weak test: 'coordinate_bnds' attribute on 'coordinate' variable is
        # required for a CF bounds variable.  This test is really testing "no
        # coordinate bounds as written by iris".
        # TODO: Fix to check variable for bounds attribute as per CF.
        expected = 0

        variables = [
            v
            for v in self.result_dataset.variables
            if "{}_bnds".format(self.coord_name)
            in self.result_dataset.variables[v].name
        ]
        actual = len(variables)

        self.assertEqual(expected, actual)

    def test_no_bounds_long_name_on_nodes(self):
        actual = self._get_variable("nodes", self.coord_name).ncattrs()

        self.assertNotIn("bounds_long_name", actual)

    def test_no_bounds_long_name_on_faces(self):
        actual = self._get_variable("face", self.coord_name).ncattrs()

        self.assertNotIn("bounds_long_name", actual)

    def test_nodes_shape(self):
        expected = (98,)

        nodes = self._get_variable("nodes", self.coord_name)
        actual = nodes.shape

        self.assertEqual(expected, actual)

    def test_faces_shape(self):
        expected = (96,)

        faces = self._get_variable("face", self.coord_name)
        actual = faces.shape

        self.assertEqual(expected, actual)

    def test_faces_type(self):
        expected = np.float64

        faces = self._get_variable("face", self.coord_name)
        actual = faces.dtype

        self.assertEqual(expected, actual)

    def test_node_type(self):
        expected = np.float64

        nodes = self._get_variable("nodes", self.coord_name)
        actual = nodes.dtype

        self.assertEqual(expected, actual)

    def test_node_units(self):
        expected = "degrees"

        nodes = self._get_variable("nodes", self.coord_name)
        actual = nodes.units

        self.assertIn(expected, actual)

    def test_face_units(self):
        expected = "degrees"

        faces = self._get_variable("face", self.coord_name)
        actual = faces.units

        self.assertIn(expected, actual)

    def test_consistent_face_node_units(self):
        faces = self._get_variable("face", self.coord_name)
        face_units = faces.units

        nodes = self._get_variable("nodes", self.coord_name)
        node_units = nodes.units

        self.assertEqual(face_units, node_units)

    def test_faces_have_correct_variable_name(self):
        name_mapping = {"longitude": "x", "latitude": "y"}
        expected = f"{TOPOLOGY_NAME}_face_{name_mapping[self.coord_name]}"

        actual = self._get_variable("face", self.coord_name).name

        self.assertEqual(actual, expected)

    def test_nodes_have_correct_variable_name(self):
        name_mapping = {"longitude": "x", "latitude": "y"}
        expected = f"{TOPOLOGY_NAME}_node_{name_mapping[self.coord_name]}"

        actual = self._get_variable("nodes", self.coord_name).name

        self.assertEqual(actual, expected)


class TestSaveLatitude(_CoordinateSaveCommon, ants.tests.TestCase):

    coord_name = "latitude"


class TestSaveLongitude(_CoordinateSaveCommon, ants.tests.TestCase):

    coord_name = "longitude"


class TestSaveMultipleCubes(ants.tests.TestCase):
    def setUp(self):
        self.cubes = ants.tests.stock.mesh_C4(load_data=True, constraint=None)

    def test_compatible_cubes(self):
        with NamedTemporaryFile("w+b") as temp_file:
            save(self.cubes, temp_file.name)

    @mock.patch("ants.fileformats._ugrid.iris_save")
    @mock.patch("ants.fileformats._ugrid.additional_save")
    def test_incompatible_cube_coordinates(self, *args):
        self.cubes[0].coord("longitude").points = (
            self.cubes[0].coord("longitude").points + 0.1
        )

        error = "Cannot save multiple UGrid cubes with different latitude or l"
        with self.assertRaisesRegex(ValueError, error):
            save(self.cubes, "foo")

    @mock.patch("ants.fileformats._ugrid.iris_save")
    @mock.patch("ants.fileformats._ugrid.additional_save")
    def test_incompatible_cube_nodes(self, *args):
        new_lats = self.cubes[0].attributes["nodes"].latitudes + 0.1
        lons = self.cubes[0].attributes["nodes"].longitudes
        self.cubes[0].attributes["nodes"] = _NodeTuple(new_lats, lons)

        error = "Cannot save cubes with different nodes attribute."
        with self.assertRaisesRegex(IOError, error):
            save(self.cubes, "foo")

    @mock.patch("ants.fileformats._ugrid.iris_save")
    @mock.patch("ants.fileformats._ugrid.additional_save")
    def test_incompatible_cube_face_node_connectivity(self, *args):
        new_data = _DataTuple(
            self.cubes[0].attributes["face_node_connectivity"].data[::-1],
            self.cubes[0].attributes["face_node_connectivity"].metadata,
        )
        self.cubes[0].attributes["face_node_connectivity"] = new_data

        error = "Cannot save cubes with different face_node_connectivity attri"
        with self.assertRaisesRegex(IOError, error):
            save(self.cubes, "foo")

    @mock.patch("ants.fileformats._ugrid.iris_save")
    @mock.patch("ants.fileformats._ugrid.additional_save")
    def test_incompatible_cube_topologies(self, *args):
        self.cubes[0].attributes["mesh_topology"].pop("face_dimension")

        error = "Cannot save cubes with different mesh_topology attri"
        with self.assertRaisesRegex(IOError, error):
            save(self.cubes, "foo")

    @mock.patch("ants.fileformats._ugrid.iris_save")
    @mock.patch("ants.fileformats._ugrid.additional_save")
    def test_incompatible_cube_names(self, *args):
        # We have special case handling for some cases with identical names,
        # but we need to catch any that slip past our name translation.
        self.cubes[0].long_name = "cube1"
        self.cubes[1].long_name = "cube1"
        cubes = iris.cube.CubeList(self.cubes)

        error = "Cannot save cubes with same name."
        with self.assertRaisesRegex(IOError, error):
            save(cubes, "foo")

    def test_data_location(self):
        expected = ("face", "face")

        with NamedTemporaryFile("w+b") as temp_file:
            save(self.cubes, temp_file.name)
            dataset = Dataset(temp_file.name)
            variable1 = dataset.variables["sample_data"]
            variable2 = dataset.variables["additional_sample_data"]
            actual = (variable1.location, variable2.location)

        self.assertEqual(expected, actual)

    def test_data_mesh(self):
        expected = (TOPOLOGY_NAME, TOPOLOGY_NAME)

        with NamedTemporaryFile("w+b") as temp_file:
            save(self.cubes, temp_file.name)
            dataset = Dataset(temp_file.name)
            variable1 = dataset.variables["sample_data"]
            variable2 = dataset.variables["additional_sample_data"]
            actual = (variable1.mesh, variable2.mesh)

        self.assertEqual(expected, actual)

    def test_data_long_names(self):
        # Requirement for XIOS is that each data variable has a long name.
        # This test goes a bit beyond that and checks the content of the long
        # names - this can be relaxed a little if needed.
        expected = ("sample_data", "additional_sample_data")

        with NamedTemporaryFile("w+b") as temp_file:
            save(self.cubes, temp_file.name)
            dataset = Dataset(temp_file.name)
            variable1 = dataset.variables["sample_data"]
            variable2 = dataset.variables["additional_sample_data"]
            actual = (variable1.long_name, variable2.long_name)

        self.assertEqual(expected, actual)

    def test_each_variable_has_correct_online_operation_attribute(self):
        # XIOS requirement, not part of CF/UGrid specs.
        expected = ("once", "once")

        with NamedTemporaryFile("w+b") as temp_file:
            save(self.cubes, temp_file.name)
            dataset = Dataset(temp_file.name)
            variable1 = dataset.variables["sample_data"]
            variable2 = dataset.variables["additional_sample_data"]
            actual = (variable1.online_operation, variable2.online_operation)

        self.assertEqual(expected, actual)
        self.assertEqual(expected, actual)


class TestSaveOtherVariables(_SaveCommon, ants.tests.TestCase):
    """Saving of variables that are neither coordinates nor data payload."""

    def _get_cf_role(self, role):
        variables = self.result_dataset.variables
        result = [
            variables[v]
            for v in variables
            if "cf_role" in variables[v].ncattrs() and variables[v].cf_role == role
        ]
        if len(result) == 0:
            self.fail("Expected cf_role {} not found".format(role))
        elif len(result) > 1:
            self.fail("Multiple cf_role attributes for {}".format(role))
        return result[0]

    def _assert_cf_role_present(self, expected):
        actual = self._get_cf_role(expected).cf_role

        self.assertEqual(expected, actual)

    def test_mesh_topology_present(self):
        self._assert_cf_role_present("mesh_topology")

    def test_topology_dimension(self):
        expected = 2

        topology = self._get_cf_role("mesh_topology")
        actual = topology.topology_dimension

        self.assertEqual(expected, actual)

    def test_topology_node_coordinates(self):
        expected = f"{TOPOLOGY_NAME}_node_y {TOPOLOGY_NAME}_node_x"

        topology = self._get_cf_role("mesh_topology")
        actual = topology.node_coordinates

        self.assertEqual(expected, actual)

    def test_topology_face_coordinates(self):
        expected = f"{TOPOLOGY_NAME}_face_y {TOPOLOGY_NAME}_face_x"

        topology = self._get_cf_role("mesh_topology")
        actual = topology.face_coordinates

        self.assertEqual(expected, actual)

    def test_topology_face_node_connectivity(self):
        expected = f"{TOPOLOGY_NAME}_face_nodes"

        topology = self._get_cf_role("mesh_topology")
        actual = topology.face_node_connectivity

        self.assertEqual(expected, actual)

    def test_topology_face_dimension(self):
        expected = "dim0"

        topology = self._get_cf_role("mesh_topology")
        actual = topology.face_dimension

        self.assertEqual(expected, actual)

    def test_topology_long_name(self):
        # From example at:
        # http://ugrid-conventions.github.io/ugrid-conventions/#3d-layered-mesh-topology
        # i.e. even for 3D layered mesh, topology name is 2D (topology is same
        # for each layer)
        expected = "Topology data of 2D unstructured mesh"

        topology = self._get_cf_role("mesh_topology")
        actual = topology.long_name

        self.assertEqual(expected, actual)

    def test_topology_dimension_name(self):
        expected = TOPOLOGY_NAME

        topology = self._get_cf_role("mesh_topology")
        actual = topology.name

        self.assertEqual(expected, actual)

    def test_face_node_mapping_present(self):
        self._assert_cf_role_present("face_node_connectivity")

    def test_face_node_mapping_datatypes(self):
        # With unittest from python 3.7 or later, can turn this into
        # parameterised test and check all mapping datatypes.  For now, settle
        # for only testing the single mapping that is required (i.e. ideally,
        # test all mappings that are present have the right dtype; but
        # settling for testing required mapping has right dtype).

        # Expected is likely to change later - large meshes known to exhaust
        # int32.  But for now, following CF spec of integers being 32 bit:
        expected = np.int32

        variables = self.result_dataset.variables
        mappings = [
            variables[v]
            for v in variables
            if (
                "cf_role" in variables[v].ncattrs()
                and "face_node_connectivity" in variables[v].cf_role
            )
        ]
        actual = mappings[0].dtype
        self.assertEqual(expected, actual)

    def test_face_node_mapping_data(self):
        expected = self.cube.attributes["face_node_connectivity"].data

        variables = self.result_dataset.variables
        mappings = [
            variables[v]
            for v in variables
            if (
                "cf_role" in variables[v].ncattrs()
                and "face_node_connectivity" in variables[v].cf_role
            )
        ]
        actual = mappings[0][:]
        self.assertArrayEqual(expected, actual)

    def test_face_node_mapping_var_name(self):
        expected = f"{TOPOLOGY_NAME}_face_nodes"
        variables = self.result_dataset.variables
        mappings = [
            variables[v]
            for v in variables
            if (
                "cf_role" in variables[v].ncattrs()
                and "face_node_connectivity" in variables[v].cf_role
            )
        ]
        self.assertEqual(1, len(mappings))
        actual = mappings[0].name
        self.assertEqual(expected, actual)

    def test_face_node_mapping_var_name_attribute_removed(self):
        expected = "var_name"

        variables = self.result_dataset.variables
        mappings = [
            variables[v]
            for v in variables
            if (
                "cf_role" in variables[v].ncattrs()
                and "face_node_connectivity" in variables[v].cf_role
            )
        ]
        self.assertEqual(1, len(mappings))
        actual = mappings[0].ncattrs()
        self.assertNotIn(expected, actual)

    def test_data_has_correct_location_attribute(self):
        expected = "face"

        actual = self.result_dataset.variables[self.data_name].location

        self.assertEqual(expected, actual)

    def test_data_has_correct_mesh_attribute(self):
        expected = TOPOLOGY_NAME

        actual = self.result_dataset.variables[self.data_name].mesh

        self.assertEqual(expected, actual)

    def test_data_has_correct_fill_value(self):
        expected = -1.7976931348623157e308

        actual = self.result_dataset.variables[self.data_name]._FillValue

        self.assertEqual(expected, actual)

    def test_fill_value_has_same_dtype_as_data(self):
        expected = self.result_dataset.variables[self.data_name].dtype

        actual = self.result_dataset.variables[self.data_name]._FillValue.dtype

        self.assertEqual(expected, actual)


class TestSaveGlobalAttributes(_SaveCommon, ants.tests.TestCase):

    # Variable attributes are tested with the variable
    def test_conventions(self):
        expected = "UGRID-1.0"

        actual = self.result_dataset.Conventions

        self.assertEqual(expected, actual)

    def test_no_face_node_connectivity(self):
        # This should only be saved as a Variable, not both a variable and an
        # attribute.
        self.assertNotIn("face_node_connectivity", self.result_dataset.ncattrs())


class TestGridMapping(_SaveCommon, ants.tests.TestCase):
    # Iris translates coordinate reference system from coords to extra
    # 'latitude_longitude' variable and 'grid_mapping' attribute on data
    # variables.  Since we don't want to handle grid mappings on nodes yet
    # (until there's a use case), we need to ensure these are not included in
    # the iris part of the save either for consistency.

    def test_no_grid_mapping_attribute(self):
        actual = self.result_dataset.variables[self.data_name].ncattrs()

        self.assertNotIn("grid_mapping", actual)

    def test_no_latitude_longitude_variable(self):
        actual = self.result_dataset.variables

        self.assertNotIn("latitude_longitude", actual)


class TestSaveMultipleDimensions(_SaveCommon, ants.tests.TestCase):
    def setUp(self):
        basis = ants.tests.stock.mesh_C4()
        self.cube = iris.cube.Cube(np.empty((2, 3, 96)))
        self.cube.add_aux_coord(basis.coord("latitude"), 2)
        self.cube.add_aux_coord(basis.coord("longitude"), 2)
        self.cube.metadata = basis.metadata
        self._temp_file = NamedTemporaryFile("w+b")
        save(self.cube, self._temp_file.name)
        self.result_dataset = Dataset(self._temp_file.name)

    def tearDown(self):
        self.result_dataset.close()
        self._temp_file.close()

    def test_face_dimension_is_last_dimension(self):
        expected = 96

        actual = self.result_dataset.dimensions["dim2"].size

        self.assertEqual(expected, actual)

    def test_face_dimension_in_topology(self):
        expected = "dim2"

        topology = self.result_dataset.variables[TOPOLOGY_NAME]
        actual = topology.face_dimension

        self.assertEqual(expected, actual)


class TestSaveFormat(_SaveCommon, ants.tests.TestCase):
    def test_netCDF_version(self):
        expected = "NETCDF4"

        actual = self.result_dataset.data_model

        self.assertTrue(actual.startswith(expected))


if __name__ == "__main__":
    ants.tests.main()
