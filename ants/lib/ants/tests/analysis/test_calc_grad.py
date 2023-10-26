# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import warnings

import ants.tests
import numpy as np
from ants.analysis import calc_grad


class Common(object):
    def test_value(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            grad_x, grad_y = calc_grad(self.source)

        # This target is simply the values as returned by the function,
        # recorded so that if behaviour results change, we are notified.
        targ_y = np.array(
            [
                [0.00000000e00, 5.77043956e-06, 1.17657014e-05],
                [0.00000000e00, 2.92269017e-06, 5.92032110e-06],
                [0.00000000e00, 7.49407735e-08, 7.49407735e-08],
            ]
        )
        targ_x = np.array(
            [
                [0.00000000e00, 0.00000000e00, 0.00000000e00],
                [-1.12411160e-06, 2.21075282e-06, -1.08664122e-06],
                [-2.99763094e-06, 5.92032110e-06, -2.92269017e-06],
            ]
        )
        self.assertTrue((np.abs(grad_x.data - targ_x) < 1e-13).all())
        self.assertTrue((np.abs(grad_y.data - targ_y) < 1e-13).all())

    def test_coordinates(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            grad_x, grad_y = calc_grad(self.source)

        reference = ants.tests.stock.geodetic(data=self.data)

        self.assertArrayEqual(
            reference.coord(axis="x").bounds, grad_x.coord(axis="x").bounds
        )
        self.assertArrayEqual(
            reference.coord(axis="x").points, grad_x.coord(axis="x").points
        )
        self.assertArrayEqual(
            reference.coord(axis="y").bounds, grad_x.coord(axis="y").bounds
        )
        self.assertArrayEqual(
            reference.coord(axis="y").points, grad_x.coord(axis="y").points
        )

        self.assertArrayEqual(
            reference.coord(axis="x").bounds, grad_y.coord(axis="x").bounds
        )
        self.assertArrayEqual(
            reference.coord(axis="x").points, grad_y.coord(axis="x").points
        )
        self.assertArrayEqual(
            reference.coord(axis="y").bounds, grad_y.coord(axis="y").bounds
        )
        self.assertArrayEqual(
            reference.coord(axis="y").points, grad_y.coord(axis="y").points
        )

        self.assertArrayEqual(
            reference.coord(axis="x").bounds, self.source.coord(axis="x").bounds
        )
        self.assertArrayEqual(
            reference.coord(axis="x").points, self.source.coord(axis="x").points
        )
        self.assertArrayEqual(
            reference.coord(axis="y").bounds, self.source.coord(axis="y").bounds
        )
        self.assertArrayEqual(
            reference.coord(axis="y").points, self.source.coord(axis="y").points
        )


class TestGlobal(Common, ants.tests.TestCase):
    def setUp(self):
        self.data = np.array([[1, 1, 1], [1, 30, 60], [1, 40, 80]])
        self.source = ants.tests.stock.geodetic(data=self.data)


class TestRotatedPole(Common, ants.tests.TestCase):
    def setUp(self):
        self.data = np.array([[1, 1, 1], [1, 30, 60], [1, 40, 80]])
        self.source = ants.tests.stock.geodetic(
            north_pole_lat=39.25, north_pole_lon=198.0, data=self.data
        )


if __name__ == "__main__":
    ants.tests.main()
