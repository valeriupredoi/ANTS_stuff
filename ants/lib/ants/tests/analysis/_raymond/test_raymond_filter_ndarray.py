# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import unittest.mock as mock

import ants.tests
import numpy as np
from ants.analysis._raymond import raymond_filter_ndarray


class TestValue(ants.tests.TestCase):
    def setUp(self):
        self.epsilon = 1.0

        length = 50
        shift = 0
        self.xvals = np.arange(shift, length + shift)
        self.perturbations = self.wave(25, length)
        self.bare = self.wave(3, length)
        self.data = self.bare + self.perturbations

    @staticmethod
    def wave(number, length):
        # Note: endpoint False as bc's are written with u_0 = u_(n+1), not
        # u_0 = u_n
        return 100 * np.cos(np.linspace(0, number * 2 * np.pi, length, endpoint=False))

    def test_1d_filter_function(self):
        # Points near the boundary exhibit a much greater error so we just
        # check that those internal ones are reasonable and that the error is
        # greater at the boundary.
        raymond_filter_ndarray(self.data, self.epsilon, self.epsilon, False)
        self.assertTrue((abs(self.data - self.bare)[7:-7] < 1).all())
        self.assertTrue((abs(self.data - self.bare) > 1).any())

    def test_1d_filter_function_periodic(self):
        raymond_filter_ndarray(self.data, self.epsilon, self.epsilon, True)
        self.assertTrue((abs(self.data - self.bare) < 1e-2).all())

    def test_2d_filter_function_periodic_x(self):
        # Ensure that when periodic is True that x-filter is periodic.
        data = np.repeat(self.data[None], 6, 0)
        bare = np.repeat(self.bare[None], 6, 0)

        raymond_filter_ndarray(data, self.epsilon, self.epsilon, True)
        self.assertTrue((abs(data - bare) < 1e-2).all())

    def test_2d_filter_function_non_periodic_x(self):
        # Ensure that when periodic is False that x-filter is not periodic.
        data = np.repeat(self.data[None], 7, 0)
        bare = np.repeat(self.bare[None], 7, 0)

        raymond_filter_ndarray(data, self.epsilon, self.epsilon, False)
        # Points near the boundary exhibit a much greater error.
        self.assertTrue((abs(data - bare)[:, 7:-7] < 1).all())
        self.assertTrue((abs(data - bare) > 1).any())

    def test_2d_filter_function_periodic_insensisitive_y(self):
        # Ensure that 'y' filtering is insensitive to being a periodic filter.
        data = np.repeat(self.data[:, None], 7, 1)
        bare = np.repeat(self.bare[:, None], 7, 1)

        for periodic in [False, True]:
            raymond_filter_ndarray(data, self.epsilon, self.epsilon, periodic)
            # Points near the boundary exhibit a much greater error.
            self.assertTrue((abs(data - bare)[:, 7:-7] < 1).all())
            self.assertTrue((abs(data - bare) > 1).any())

    def test_non_periodic_varying_epsilon(self):
        # Test the isotropic case.
        # We pass a list of epsilon all of the same value to check that it
        # returns the same results as the anisotropic case (constant epsilon).
        tar_data = self.data.copy()
        raymond_filter_ndarray(tar_data, self.epsilon, self.epsilon, periodic=False)

        eps_row = np.empty(self.data.shape[0])
        eps_row.fill(self.epsilon)
        raymond_filter_ndarray(self.data, self.epsilon, eps_row, periodic=False)
        self.assertArrayAlmostEqual(self.data, tar_data)

    def test_periodic_varying_epsilon(self):
        # Test the isotropic case.
        # We pass a list of epsilon all of the same value to check that it
        # returns the same results as the anisotropic case (constant epsilon).
        tar_data = self.data.copy()
        raymond_filter_ndarray(tar_data, self.epsilon, self.epsilon, periodic=True)

        eps_row = np.empty(self.data.shape[0])
        eps_row.fill(self.epsilon)
        raymond_filter_ndarray(self.data, self.epsilon, eps_row, periodic=True)

        self.assertArrayAlmostEqual(self.data, tar_data)


class TestCall(ants.tests.TestCase):
    def setUp(self):
        patch = mock.patch("ants.analysis._raymond.filter_rows_constant_epsilon")
        self.patch_ceps = patch.start()
        self.addCleanup(patch.stop)
        patch = mock.patch("ants.analysis._raymond.filter_rows_varying_epsilon")
        self.patch_veps = patch.start()
        self.addCleanup(patch.stop)
        patch = mock.patch("ants.analysis._raymond.filter_columns")
        self.patch_cols = patch.start()
        self.data = mock.Mock(shape=(3, 4), ndim=2)

    def test_varying_epsilon_call_args(self):
        eps_row = np.ones(self.data.shape[0])
        raymond_filter_ndarray(
            self.data,
            mock.sentinel.epsilon_col,
            eps_row,
            periodic=mock.sentinel.periodic,
        )
        self.assertFalse(self.patch_ceps.called)
        self.patch_veps.assert_called_once_with(
            self.data, eps_row, mock.sentinel.periodic
        )
        self.patch_cols.assert_called_once_with(self.data, mock.sentinel.epsilon_col)

    def test_constant_epsilon_call_args(self):
        eps_row = 1
        raymond_filter_ndarray(
            self.data,
            mock.sentinel.epsilon_col,
            eps_row,
            periodic=mock.sentinel.periodic,
        )
        self.patch_ceps.assert_called_once_with(
            self.data, eps_row, mock.sentinel.periodic
        )
        self.assertFalse(self.patch_veps.called)
        self.patch_cols.assert_called_once_with(self.data, mock.sentinel.epsilon_col)


if __name__ == "__main__":
    ants.tests.main()
