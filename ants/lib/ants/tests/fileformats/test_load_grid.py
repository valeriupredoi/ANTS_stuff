# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import tempfile
import unittest.mock as mock

import ants
import ants.tests
import iris
import numpy as np
from ants.fileformats import load_grid


@mock.patch("ants.fileformats._ugrid.load_mesh", side_effect=RuntimeError)
class TestAll(ants.tests.TestCase):
    def setUp(self):
        patch = mock.patch("ants.fileformats.load")
        self.mock_ants_load = patch.start()
        self.addCleanup(patch.stop)

        self.surface_altitude = ants.tests.stock.geodetic((3, 4))
        self.surface_altitude.data = self.surface_altitude.data + 1000
        self.surface_altitude.standard_name = "surface_altitude"
        self.surface_altitude.units = "m"

        self.model_vertical = iris.cube.Cube(np.zeros(3))
        mln = iris.coords.DimCoord(
            np.arange(1, 4),
            standard_name="model_level_number",
            units="1",
            attributes={"positive": "up"},
        )
        sigma = iris.coords.AuxCoord(
            [0.2, 0.4, 0.6],
            bounds=[[0.1, 0.3], [0.3, 0.5], [0.5, 0.7]],
            long_name="sigma",
            units="1",
        )
        level_height = iris.coords.AuxCoord(
            [0.5, 1.0, 1.5],
            bounds=[[0.25, 0.75], [0.75, 1.25], [1.25, 1.75]],
            long_name="level_height",
            attributes={"positive": "up"},
            units="m",
        )
        self.model_vertical.add_dim_coord(mln, 0)
        self.model_vertical.add_aux_coord(sigma, 0)
        self.model_vertical.add_aux_coord(level_height, 0)
        self.model_vertical.long_name = "Model vertical definition"

        self.model_horizontal = ants.tests.stock.geodetic((3, 4))
        self.model_horizontal.long_name = "Model Grid"

    def assert_cube_equal(self, result, target):
        # Replace data as lazy
        result.data = target.data
        target.rename(None)
        target.units = None
        self.assertEqual(result, target)

    def test_horizontal(self, _):
        self.mock_ants_load.return_value = iris.cube.CubeList(
            [self.model_horizontal.copy()]
        )
        res = load_grid("dummy_filename")
        self.assert_cube_equal(res, self.model_horizontal)

    def test_vertical(self, _):
        self.mock_ants_load.return_value = iris.cube.CubeList(
            [self.model_vertical.copy()]
        )
        res = load_grid("dummy_filename")
        self.assert_cube_equal(res, self.model_vertical)

    def test_surface_altitude(self, _):
        self.mock_ants_load.return_value = iris.cube.CubeList(
            [self.surface_altitude.copy()]
        )
        res = load_grid("dummy_filename")
        self.assertArrayEqual(
            res.coord("surface_altitude").points, self.surface_altitude.data
        )
        self.assertCML(
            res, ("fileformats", "load_grid", "surface_altitude.cml"), checksum=False
        )

    def test_altitude(self, _):
        # Hybrid height factory check.
        self.mock_ants_load.return_value = iris.cube.CubeList(
            [self.model_vertical.copy(), self.surface_altitude.copy()]
        )
        res = load_grid("dummy_filename")
        self.assertCML(
            res, ("fileformats", "load_grid", "altitude.cml"), checksum=False
        )

    def test_altitude_alt_mapping(self, _):
        self.surface_altitude.transpose((1, 0))
        self.mock_ants_load.return_value = iris.cube.CubeList(
            [self.model_vertical.copy(), self.surface_altitude.copy()]
        )
        res = load_grid("dummy_filename")
        self.assertCML(
            res,
            ("fileformats", "load_grid", "altitude_alt_mapping.cml"),
            checksum=False,
        )

    def test_general_target_field_no_vertical(self, _):
        # Ensure that we can extract horizontal information from any random
        # cube picked up with any mapping (not just surface altitude or our
        # namelist grids).
        # Also ensure that we inherrit the metadata in this case.
        horizontal = ants.tests.stock.geodetic((5, 3, 4))
        horizontal.transpose((1, 0, 2))
        self.mock_ants_load.return_value = iris.cube.CubeList([horizontal])
        res = load_grid("dummy_filename")
        self.assert_cube_equal(res, horizontal)

    def test_general_target_field_with_vertical(self, _):
        # Ensure that we can extract horizontal information from any random
        # cube picked up with any mapping and associate vertical components to
        # it.
        horizontal = ants.tests.stock.geodetic((5, 3, 4))
        horizontal.transpose((1, 0, 2))
        self.mock_ants_load.return_value = iris.cube.CubeList(
            [horizontal, self.model_vertical.copy()]
        )
        res = load_grid("dummy_filename")
        self.assertCML(
            res, ("fileformats", "load_grid", "horizontal_vertical.cml"), checksum=False
        )

    def test_remove_disk_dependence(self, _):
        self.mock_ants_load.return_value = iris.cube.CubeList(
            [self.model_horizontal.copy()]
        )
        res = load_grid("dummy_filename")
        self.assertTrue(res.has_lazy_data())

    def test_multiple_identical_grids(self, _):
        self.mock_ants_load.return_value = iris.cube.CubeList(
            [self.model_horizontal.copy(), self.model_horizontal.copy()]
        )
        res = load_grid("dummy_filename")
        self.assert_cube_equal(res, self.model_horizontal)

    def test_multiple_conflicting_grids(self, _):
        horizontal_grid_alt = self.model_horizontal.copy()
        horizontal_grid_alt.coord(axis="x").points = (
            horizontal_grid_alt.coord(axis="x").points + 1
        )
        self.mock_ants_load.return_value = iris.cube.CubeList(
            [self.model_horizontal.copy(), horizontal_grid_alt]
        )

        msg = "Unable to resolve a single grid, conflicting coordinate: " "longitude"
        with self.assertRaisesRegex(ValueError, msg):
            load_grid("dummy_filename")

    def test_no_grid(self, _):
        self.mock_ants_load.return_value = None
        msg = "no cubes found"
        with self.assertRaisesRegex(iris.exceptions.ConstraintMismatchError, msg):
            load_grid("dummy_filename")


@mock.patch("ants.fileformats._ugrid.load_mesh", side_effect=RuntimeError)
class TestVariableResolution(ants.tests.TestCase):
    def setUp(self):
        regular_namelist = b"""&GRID
POINTS_LAMBDA_TARG=2,POINTS_PHI_TARG=2,
PHI_POLE=37.5,LAMBDA_POLE=177.5,ROTATED=T
/
"""
        variable_namelist = b"""&HORIZGRID
 LAMBDA_INPUT_P=1, 2, 3
 LAMBDA_INPUT_U=1.5, 2.5, 3.5
 PHI_INPUT_P=10, 20, 30
 PHI_INPUT_V=5, 15, 25
/
"""
        with tempfile.NamedTemporaryFile() as regular_fh:
            with tempfile.NamedTemporaryFile() as variable_fh:
                regular_fh.write(regular_namelist)
                variable_fh.write(variable_namelist)
                regular_fh.seek(0)
                variable_fh.seek(0)
                self.result = load_grid([regular_fh.name, variable_fh.name])

        ellipsoid = iris.coord_systems.GeogCS(6371229.0)
        self.expected_coord_system = iris.coord_systems.RotatedGeogCS(
            37.5, 177.5, ellipsoid=ellipsoid
        )

    def test_returns_cube(self, _):
        self.assertIsInstance(self.result, iris.cube.Cube)

    def test_latitude(self, _):
        expected = iris.coords.DimCoord(
            [10, 20, 30],
            bounds=[[5, 15], [15, 25], [25, 35]],
            units="degrees",
            standard_name="grid_latitude",
            coord_system=self.expected_coord_system,
        )

        actual = self.result.coord(axis="y")

        self.assertEqual(actual, expected)

    def test_longitude(self, _):
        expected = iris.coords.DimCoord(
            [1, 2, 3],
            bounds=[[0.5, 1.5], [1.5, 2.5], [2.5, 3.5]],
            units="degrees",
            standard_name="grid_longitude",
            coord_system=self.expected_coord_system,
        )

        actual = self.result.coord(axis="x")

        self.assertEqual(actual, expected)


class TestMesh(ants.tests.TestCase):
    # Mesh tests can't be in TestAll - the mock in setUp() breaks mesh loading
    # (needed for ants.tests.stock.mesh_C4):
    def test_mesh_only(self):
        mesh = ants.tests.stock.mesh_C4()
        expected = (mesh.coord(axis="x"), mesh.coord(axis="y"))

        mesh_path = ants.tests.get_data_path(("stock", "mesh_C4.nc"))
        result = load_grid(mesh_path)
        actual = (result.coord(axis="x"), result.coord(axis="y"))

        self.assertEqual(actual, expected)

    def test_mesh_with_levels(self):
        template = ants.tests.stock.simple_4d_with_hybrid_height()
        vertical = iris.cube.Cube(np.zeros(template.shape[1]))
        vertical.add_dim_coord(template.coord("model_level_number"), 0)
        vertical.add_aux_coord(template.coord("level_height"), 0)
        vertical.add_aux_coord(template.coord("sigma"), 0)

        mesh = ants.tests.stock.mesh_C4()
        expected = (
            mesh.coord(axis="x"),
            mesh.coord(axis="y"),
            vertical.coord(axis="z"),
        )

        vertical = iris.cube.CubeList((vertical,))

        mesh_path = ants.tests.get_data_path(("stock", "mesh_C4.nc"))
        with mock.patch("ants.fileformats.load", return_value=vertical):
            # Shortcircuit error thrown in load_mesh - avoids multiple
            # additional mocks:
            with mock.patch(
                "ants.fileformats._ugrid.iris.load",
                side_effect=(
                    (
                        iris.load(mesh_path),  # UGrid file first...
                        RuntimeError(),  # Then error for regular file
                    )
                ),
            ):
                result = load_grid((mesh_path, "foo"))

        actual = (
            result.coord(axis="x"),
            result.coord(axis="y"),
            result.coord(axis="z"),
        )

        self.assertEqual(actual, expected)

    def test_mesh_attributes(self):
        mesh = ants.tests.stock.mesh_C4()
        expected = mesh.attributes

        mesh_path = ants.tests.get_data_path(("stock", "mesh_C4.nc"))
        result = load_grid(mesh_path)
        actual = result.attributes

        self.assertEqual(actual, expected)


if __name__ == "__main__":
    ants.tests.main()
