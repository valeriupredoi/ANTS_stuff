# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
import iris
import numpy as np
from ants.regrid.rectilinear import _fill_outside_bounds


class Common(object):
    def setUp(self):
        self.source = ants.tests.stock.geodetic((5, 5), ylim=(-90, 30), xlim=(-160, 80))
        self.target = ants.tests.stock.geodetic((5, 5))
        self.target.data = self.target.data.astype("float32")

        self.result = np.array(
            [
                [np.nan, 1, 2, np.nan, np.nan],
                [np.nan, 6, 7, np.nan, np.nan],
                [np.nan, 11, 12, np.nan, np.nan],
                [np.nan, np.nan, np.nan, np.nan, np.nan],
                [np.nan, np.nan, np.nan, np.nan, np.nan],
            ],
            dtype=np.float32,
        )


class Test2D(Common, ants.tests.TestCase):
    def test_increasing(self):
        _fill_outside_bounds(self.source, self.target, np.NaN)
        self.assertArrayEqual(self.target.data, self.result)

    @staticmethod
    def _invert_coord(coord):
        coord.points = coord.points[::-1]
        coord.bounds = coord.bounds[::-1, ::-1]

    def test_decreasing_source(self):
        sx = self.source.coord(axis="x")
        self._invert_coord(sx)
        sy = self.source.coord(axis="y")
        self._invert_coord(sy)

        _fill_outside_bounds(self.source, self.target, np.NaN)
        self.assertArrayEqual(self.target.data, self.result)

    def test_decreasing_target(self):
        # We copy the data here due to a numpy bug
        # https://github.com/numpy/numpy/issues/8264
        self.target.data = self.target.data[::-1, ::-1].copy()

        tx = self.target.coord(axis="x")
        self._invert_coord(tx)
        ty = self.target.coord(axis="y")
        self._invert_coord(ty)

        _fill_outside_bounds(self.source, self.target, np.NaN)
        self.assertArrayEqual(self.target.data, self.result[::-1, ::-1])

    def test_masked(self):
        # Ensure that masked elements outside the source extent become nan and
        # unmasked while masked elements inside the extent remain masked.
        self.target.data = np.ma.array(self.target.data)
        self.target.data[0, :] = np.ma.masked
        _fill_outside_bounds(self.source, self.target, np.NaN)

        self.result = np.ma.array(self.result)
        self.result[0, 1:3] = np.ma.masked
        self.assertMaskedArrayEqual(self.target.data, self.result)


class TestND(Common, ants.tests.TestCase):
    def test_zyx_source(self):
        # Ensure no sensitivity to source dimension mapping
        cube1 = self.source
        cube2 = self.source.copy()

        coord = iris.coords.DimCoord(0, long_name="bla")
        cube1.add_aux_coord(coord, None)
        coord = iris.coords.DimCoord(1, long_name="bla")
        cube2.add_aux_coord(coord, None)

        self.source = iris.cube.CubeList([cube1, cube2]).merge_cube()

        _fill_outside_bounds(self.source, self.target, np.NaN)
        self.assertArrayEqual(self.target.data, self.result)

    def test_xzy(self):
        # Ensure no sensitivity to target dimension mapping
        cube1 = self.target
        cube2 = self.target.copy()

        coord = iris.coords.DimCoord(0, long_name="bla")
        cube1.add_aux_coord(coord, None)
        coord = iris.coords.DimCoord(1, long_name="bla")
        cube2.add_aux_coord(coord, None)

        self.target = iris.cube.CubeList([cube1, cube2]).merge_cube()

        self.result = np.vstack([self.result[None, ...], self.result[None, ...]])
        self.target.transpose((2, 0, 1))
        self.result = self.result.transpose((2, 0, 1))

        _fill_outside_bounds(self.source, self.target, np.NaN)
        self.assertArrayEqual(self.target.data, self.result)


if __name__ == "__main__":
    ants.tests.main()
