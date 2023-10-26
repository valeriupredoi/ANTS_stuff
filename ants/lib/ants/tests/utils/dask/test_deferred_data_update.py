# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import unittest.mock as mock

import ants.tests
import numpy as np
from ants.utils.dask import deferred_data_update


class TestValues(ants.tests.TestCase):
    def test_2d_newdata(self):
        arr = np.arange(16).reshape((4, 4))
        tarr = np.array([[100, 200], [300, 400]])
        slices = tuple([slice(0, 2), slice(1, 3)])
        result = deferred_data_update(arr, tarr, slices)

        tar = np.array(
            [[0, 100, 200, 3], [4, 300, 400, 7], [8, 9, 10, 11], [12, 13, 14, 15]]
        )
        self.assertArrayEqual(result, tar)


class TestExceptions(ants.tests.TestCase):
    def test_nd_source(self):
        source = mock.Mock(name="data", ndim=1)
        target = mock.Mock(name="data", ndim=2)
        msg = "Expected 2D source data, got 1 instead"
        with self.assertRaisesRegex(ValueError, msg):
            deferred_data_update(source, target, mock.sentinel.slices)

    def test_nd_target(self):
        source = mock.Mock(name="data", ndim=2)
        target = mock.Mock(name="data", ndim=1)
        msg = "Expected 2D target data, got 1 instead"
        with self.assertRaisesRegex(ValueError, msg):
            deferred_data_update(source, target, mock.sentinel.slices)


if __name__ == "__main__":
    ants.tests.main()
