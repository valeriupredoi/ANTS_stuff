# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
from unittest.mock import patch

import ants.tests
import iris
import numpy as np
from ants.regrid._ugrid import _UGridRegridder as UGridRegridder


@ants.tests.skip_esmf
@patch("ants.regrid._ugrid._build_mesh")
@patch("ants.regrid._ugrid._build_grid")
class TestInit(ants.tests.TestCase):
    def test_build_grid_call(self, build_grid, _):
        UGridRegridder(source="foo", target="bar")
        build_grid.assert_called_once_with("foo")

    def test_build_mesh_call(self, _, build_mesh):
        UGridRegridder(source="foo", target="bar")
        build_mesh.assert_called_once_with("bar")


@ants.tests.skip_esmf
class TestCall(ants.tests.TestCase):
    def setUp(self):
        self.mesh = ants.tests.stock.mesh_C4()
        self.regrid = None

    def tearDown(self):
        # Need to wipe cache to preserve test isolation
        self.regrid.manager.clear()

    @patch("ants.regrid._ugrid._create_result")
    def test_create_result_call(self, create_result):
        grid = ants.tests.stock.geodetic((3, 3))
        self.regrid = UGridRegridder(source=grid, target=self.mesh)
        self.regrid(ants.tests.stock.geodetic((3, 3)))
        create_result.assert_called_once_with(grid, self.mesh)

    def test_reject_cube_with_unsupported_dimension_order(self):
        grid = ants.tests.stock.simple_3d_time_varying()
        grid.transpose((2, 1, 0))
        self.regrid = UGridRegridder(source=grid, target=self.mesh)

        with self.assertRaisesRegex(RuntimeError, "Unsupported dimension"):
            self.regrid(grid)

    @patch("ants.regrid._ugrid.np.isnan", return_value=False)
    @patch("ants.regrid._ugrid._create_result")
    @patch("ants.regrid._ugrid.CachedRegridManager._create_regrid")
    def test_regrid_cached(self, mock_regrid, *args):
        grid = ants.tests.stock.geodetic((3, 3))
        self.regrid = UGridRegridder(source=grid, target=self.mesh)
        self.regrid(grid)
        grid = ants.tests.stock.geodetic((3, 3))
        mesh = ants.tests.stock.mesh_C4()
        UGridRegridder(source=grid, target=mesh)(grid)
        mock_regrid.assert_called_once()

    @patch("ants.regrid._ugrid.np.isnan", return_value=False)
    @patch("ants.regrid._ugrid._create_result")
    @patch("ants.regrid._ugrid.CachedRegridManager._create_regrid")
    def test_regrid_uncached(self, mock_regrid, *args):
        grid1 = ants.tests.stock.geodetic((3, 3))
        grid2 = ants.tests.stock.geodetic((4, 3))
        self.regrid = UGridRegridder(source=grid1, target=self.mesh)
        self.regrid(grid1)
        UGridRegridder(source=grid2, target=self.mesh)(grid2)
        self.assertEqual(2, mock_regrid.call_count)


@ants.tests.skip_esmf
class TestMultiples(ants.tests.TestCase):
    def setUp(self):
        self.mesh = ants.tests.stock.mesh_C4()
        self.regrid = None

    def test_consistent_result_for_multiple_uncached_calls(self):
        # Values tuned to be as low as possible, while still showing test
        # failure reasonably consistently.
        N = 10
        grid = ants.tests.stock.geodetic((10, 10))
        ref = None
        for i in range(N):
            self.regrid = UGridRegridder(source=grid, target=self.mesh)
            result = self.regrid(grid).data
            if ref is None:
                ref = result
            self.assertArrayEqual(result, ref)
            # Clears cache
            self.regrid.manager.clear()


@ants.tests.skip_esmf
class TestNaNBehaviour(ants.tests.TestCase):
    def setUp(self):
        self.mesh = ants.tests.stock.mesh_C4()
        self.regrid = None
        self.mock_result = np.ma.ones(self.mesh.data.shape, dtype=np.float)
        self.mock_result[0] = np.nan
        self.grid = ants.tests.stock.geodetic((3, 3))
        self.grid.data = self.grid.data.astype(np.float)

    def tearDown(self):
        # Need to wipe cache to preserve test isolation
        self.regrid.manager.clear()

    def assertDoesNotRaise(self, source):
        self.regrid = UGridRegridder(source=source, target=self.mesh)
        with patch(
            "ants.regrid._ugrid._UGridRegridder._regrid_2d_field",
            return_value=self.mock_result,
        ):
            result = self.regrid(source)
        self.assertIsNotNone(result)

    def test_valid_data(self):
        self.regrid = UGridRegridder(source=self.grid, target=self.mesh)
        self.assertIsNotNone(self.regrid(ants.tests.stock.geodetic((3, 3))))

    def test_source_has_NaNs(self):
        source = self.grid.copy(data=self.grid.data)
        source.data[0, 0] = np.nan
        self.assertDoesNotRaise(source)

    def test_source_has_masked_values(self):
        source = self.grid.copy(data=self.grid.data)
        source.data[0, 0] = np.ma.masked
        self.assertDoesNotRaise(source)

    def test_source_is_local(self):
        crs = iris.coord_systems.GeogCS(6371229.0)
        grid = ants.tests.stock.gen_regular_cube(shape=(3, 3), xlim=(0, 180), crs=crs)
        source = grid.copy(data=grid.data)
        self.assertDoesNotRaise(source)

    def test_valid_source_invalid_result(self):
        self.regrid = UGridRegridder(source=self.grid, target=self.mesh)
        with patch(
            "ants.regrid._ugrid._UGridRegridder._regrid_2d_field",
            return_value=self.mock_result,
        ):
            with self.assertRaisesRegex(
                RuntimeError, "Cube unknown contains unexpected NaN values."
            ):
                self.regrid(ants.tests.stock.geodetic((3, 3)))


if __name__ == "__main__":
    ants.tests.main()
