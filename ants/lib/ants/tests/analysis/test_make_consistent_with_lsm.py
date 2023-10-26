# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import unittest.mock as mock

import ants.tests
import numpy as np
from ants.analysis import make_consistent_with_lsm


class TestRegular(ants.tests.TestCase):
    def setUp(self):
        self.source = ants.tests.stock.geodetic((2, 2))
        self.lsm = ants.tests.stock.geodetic((2, 2))
        self.lsm.data = self.lsm.data.astype("bool")

    @mock.patch("ants.analysis.FillMissingPoints")
    def test_filler(self, patch_fill):
        make_consistent_with_lsm(self.source, self.lsm, False)
        patch_fill.assert_called_once_with(self.source, target_mask=self.lsm)

    @mock.patch("ants.utils.cube.guess_horizontal_bounds")
    def test_guess_bounds(self, patch_guess):
        make_consistent_with_lsm(self.source, self.lsm, False)
        patch_guess.assert_has_calls([mock.call(self.lsm), mock.call(self.source)])
        self.assertEqual(2, patch_guess.call_count)


class TestUgrid(ants.tests.TestCase):
    def setUp(self):
        self.source = ants.tests.stock.mesh_C4(
            load_data=True, constraint="sample_data"
        )[0]
        self.lsm = self.source.copy(data=np.zeros(self.source.shape, dtype="bool"))
        self.lsm.data[::3] = 1

    @mock.patch("ants.analysis._UGridFillMissingPoints")
    def test_filler(self, patch_fill):
        make_consistent_with_lsm(self.source, self.lsm, False)
        patch_fill.assert_called_once_with(self.source, target_mask=self.lsm)

    @mock.patch("ants.utils.cube.guess_horizontal_bounds")
    def test_no_guess_bounds(self, patch_guess):
        make_consistent_with_lsm(self.source, self.lsm, False)
        patch_guess.assert_not_called()


class TestLSMWarning(ants.tests.TestCase):
    def setUp(self):
        self.source = ants.tests.stock.geodetic((2, 2))
        self.lsm = ants.tests.stock.geodetic((2, 2))
        self.lsm.data = self.lsm.data.astype("bool")

    def test_no_warning(self):
        with mock.patch("warnings.warn") as mock_warn:
            make_consistent_with_lsm(self.source, self.lsm, False)
        self.assertEqual(mock_warn.call_count, 0)

    def test_masked_value(self):
        self.lsm.data = np.ma.masked_array(self.lsm.data, mask=[[0, 1], [0, 0]])
        with mock.patch("warnings.warn") as mock_warn:
            make_consistent_with_lsm(self.source, self.lsm, False)
        self.assertEqual(mock_warn.call_count, 1)

    def test_non_bool_value(self):
        self.lsm.data = np.array([[0, 1], [2, 0]])
        with mock.patch("warnings.warn") as mock_warn:
            make_consistent_with_lsm(self.source, self.lsm, False)
        self.assertEqual(mock_warn.call_count, 1)

    def test_non_bool_value_between_0_and_1(self):
        self.lsm.data = np.array([[0, 1], [0.5, 0]])
        with mock.patch("warnings.warn") as mock_warn:
            make_consistent_with_lsm(self.source, self.lsm, False)
        self.assertEqual(mock_warn.call_count, 1)

    def test_masked_value_and_non_bool_value(self):
        self.lsm.data = np.array([[0, 1], [2, 0]])
        self.lsm.data = np.ma.masked_array(self.lsm.data, mask=[[0, 1], [0, 0]])
        with mock.patch("warnings.warn") as mock_warn:
            make_consistent_with_lsm(self.source, self.lsm, False)
        self.assertEqual(mock_warn.call_count, 2)


if __name__ == "__main__":
    ants.tests.main()
