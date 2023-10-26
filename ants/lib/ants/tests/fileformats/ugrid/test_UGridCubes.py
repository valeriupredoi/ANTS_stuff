# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
import iris
import numpy as np
from ants.fileformats._ugrid import _UGridCubes as UGridCubes


class _UGridCubesCommon(object):
    def _get_raw_cube(self, var_name):
        result = [cube for cube in self.raw_cubes if var_name == cube.var_name]
        assert len(result) == 1
        return result[0]

    def _get_coordinate_cube(self, attribute_name, standard_name):
        coordinate = [
            attribute
            for attribute in self.expected_topology.attributes[attribute_name].split()
            if self._get_raw_cube(attribute).standard_name == standard_name
        ]
        assert len(coordinate) == 1
        results = [cube for cube in self.raw_cubes if cube.var_name == coordinate[0]]
        assert len(results) == 1
        return results[0]

    def test_node_longitudes(self):
        expected = self._get_coordinate_cube("node_coordinates", "longitude").data

        actual = self.ugrid_cubes._construct_mesh().attributes.get("nodes").longitudes

        self.assertArrayEqual(expected, actual)

    def test_node_latitudes(self):
        expected = self._get_coordinate_cube("node_coordinates", "latitude").data

        actual = self.ugrid_cubes._construct_mesh().attributes.get("nodes").latitudes

        self.assertArrayEqual(expected, actual)

    def test_error_for_multiple_topologies(self):
        topology = [
            c for c in self.raw_cubes if c.attributes.get("cf_role") == "mesh_topology"
        ][0]
        self.raw_cubes.append(topology)
        ugrid = UGridCubes(self.raw_cubes)
        with self.assertRaisesRegex(
            RuntimeError, "Expected exactly 1 topology, found 2."
        ):
            ugrid._get_topologies()


class TestMesh(_UGridCubesCommon, ants.tests.TestCase):
    def setUp(self):
        self.raw_cubes = iris.load(ants.tests.get_data_path(("stock", "mesh_C4.nc")))
        self.ugrid_cubes = UGridCubes(self.raw_cubes)
        expected_topology = [
            c for c in self.raw_cubes if c.attributes.get("cf_role") == "mesh_topology"
        ]
        assert len(expected_topology) == 1
        self.expected_topology = expected_topology[0]

    def test_face_longitudes(self):
        expected = self._get_coordinate_cube("face_coordinates", "longitude").data

        actual = self.ugrid_cubes._construct_mesh().coord("longitude").points

        self.assertArrayEqual(expected, actual)

    def test_face_latitudes(self):
        expected = self._get_coordinate_cube("face_coordinates", "latitude").data

        actual = self.ugrid_cubes._construct_mesh().coord("latitude").points

        self.assertArrayEqual(expected, actual)


class TestData(_UGridCubesCommon, ants.tests.TestCase):
    def setUp(self):
        self.raw_cubes = iris.load(ants.tests.get_data_path(("stock", "data_C4.nc")))
        self.ugrid_cubes = UGridCubes(self.raw_cubes, data_constraint="sample_data")
        expected_topology = [
            c for c in self.raw_cubes if c.attributes.get("cf_role") == "mesh_topology"
        ]
        assert len(expected_topology) == 1
        self.expected_topology = expected_topology[0]

    def test_face_longitudes(self):
        expected = (
            self.raw_cubes.extract("sample_data", strict=True).coord("longitude").points
        )
        actual = self.ugrid_cubes._construct_mesh().coord("longitude").points

        self.assertArrayEqual(expected, actual)

    def test_face_latitudes(self):
        expected = (
            self.raw_cubes.extract("sample_data", strict=True).coord("latitude").points
        )

        actual = self.ugrid_cubes._construct_mesh().coord("latitude").points

        self.assertArrayEqual(expected, actual)


class TestConstructCoordinate(ants.tests.TestCase):
    def setUp(self):
        """
        Setup of the initial test data.

        From the following, understanding the expected values in each test
        should be trivial.

        4 --- 3 --- 6

        |     |     |
        | 101 | 102 |
        |     |     |

        1 --- 2 --- 5

        Arbitrarily defining node latitudes as 1 for nodes 1, 2 and 5, and 3
        for nodes 4, 3 and 6.

        Defining node longitudes as 1 for nodes 1 and 4, 3 for nodes 2 and 3,
        and 5 for nodes 5 and 6.

        Face latitudes are 2 for all faces.

        Face longitudes are 2 for face 101, and 4 for face 102.

        Face to node mapping is [1, 2, 3, 4] for face 101, and [2, 5, 6, 3]
        for face 102 (note: UGRID spec mandates counterclockwise mapping).

        Real meshes use 1-based - rather than 0-based - indexing for the
        face_node_mapping.  It's imperative that this test data preserves
        that.  Everything else in this test data is arbitrary.

        """
        self.node_lats = iris.cube.Cube(
            (1, 1, 3, 3, 1, 3),
            long_name="latitude of 2D mesh nodes.",
        )
        self.node_lons = iris.cube.Cube(
            (1, 3, 3, 1, 5, 5),
            long_name="longitude of 2D mesh nodes.",
        )
        self.face_lats = iris.cube.Cube(np.array((2, 2)))
        self.face_lons = iris.cube.Cube(np.array((2, 4)))

        _mapping = "Maps every quadrilateral face to its four corner nodes."
        self.face_node_mapping = iris.cube.Cube(
            np.array([[1, 2, 3, 4], [2, 5, 6, 3]]),
            long_name=_mapping,
        )
        self.raw_cubes = iris.load(ants.tests.get_data_path(("stock", "data_C4.nc")))
        self.ugrid_cubes = UGridCubes(self.raw_cubes)

    def test_CF_name(self):
        expected = "latitude"

        coord = self.ugrid_cubes._construct_coordinate(
            self.node_lats, self.face_lats, self.face_node_mapping, "latitude"
        )
        actual = coord.standard_name

        self.assertEqual(actual, expected)

    def test_latitude_type(self):
        expected = iris.coords.AuxCoord

        actual = self.ugrid_cubes._construct_coordinate(
            self.node_lats, self.face_lats, self.face_node_mapping
        )

        self.assertIsInstance(actual, expected)

    def test_latitude_points(self):
        expected = np.array((2, 2))

        actual = self.ugrid_cubes._construct_coordinate(
            self.node_lats, self.face_lats, self.face_node_mapping
        ).points

        self.assertArrayEqual(actual, expected)

    def test_latitude_bounds(self):
        expected = np.array([[1, 1, 3, 3], [1, 1, 3, 3]])

        actual = self.ugrid_cubes._construct_coordinate(
            self.node_lats, self.face_lats, self.face_node_mapping
        ).bounds

        self.assertArrayEqual(actual, expected)

    def test_longitude_type(self):
        expected = iris.coords.AuxCoord

        actual = self.ugrid_cubes._construct_coordinate(
            self.node_lons, self.face_lons, self.face_node_mapping
        )

        self.assertIsInstance(actual, expected)

    def test_longitude_points(self):
        expected = np.array((2, 4))

        actual = self.ugrid_cubes._construct_coordinate(
            self.node_lons, self.face_lons, self.face_node_mapping
        ).points

        self.assertArrayEqual(actual, expected)

    def test_longitude_bounds(self):
        expected = np.array([[1, 3, 3, 1], [3, 5, 5, 3]])

        actual = self.ugrid_cubes._construct_coordinate(
            self.node_lons, self.face_lons, self.face_node_mapping
        ).bounds

        self.assertArrayEqual(actual, expected)


if __name__ == "__main__":
    ants.tests.main()
