# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import copy
import unittest.mock as mock

import ants.tests
from ants.analysis.filters import raymond


class TestAll(ants.tests.TestCase):
    def setUp(self):
        patch = mock.patch("ants.analysis._raymond.raymond_filter_ndarray")
        self.patch_raymond_filt_ndarray = patch.start()
        self.addCleanup(patch.stop)

        self.source = ants.tests.stock.geodetic((3, 4))
        self.source.coord(axis="x").circular = False

    def assert_raymond_called_with(self, *args):
        args = list(args)
        for i in range(len(args)):
            arg = args[i]
            if hasattr(arg, "__iter__"):
                self.assertArrayAlmostEqual(
                    self.patch_raymond_filt_ndarray.call_args_list[0][0][i], arg
                )
            else:
                self.assertEqual(
                    self.patch_raymond_filt_ndarray.call_args_list[0][0][i], arg
                )

    def test_filter_length_scale(self):
        raymond(self.source, filter_length_scale=1, isotropic=False)
        self.assert_raymond_called_with(
            self.source.data, 449.3211491564997, 449.3211491564997, False
        )

    def test_filter_length_scale_isotropic(self):
        raymond(self.source, filter_length_scale=1, isotropic=True)
        self.assert_raymond_called_with(
            self.source.data,
            449.3211491564997,
            [0.03794428, 0.03528191, 0.03794428],
            False,
        )

    def test_epsilon(self):
        raymond(self.source, epsilon=1, isotropic=False)
        self.assert_raymond_called_with(self.source.data, 1, 1, False)

    def test_epsilon_isotropic(self):
        raymond(self.source, epsilon=1, isotropic=True)
        self.assert_raymond_called_with(
            self.source.data, 1, [11.236847, 1.0, 11.236847], False
        )

    def test_circular(self):
        # Check periodic.
        self.source.coord(axis="x").circular = True
        raymond(self.source, epsilon=1, isotropic=False)
        self.assert_raymond_called_with(self.source.data, 1, 1, True)

    def test_isotropic_not_supported(self):
        self.source.coord(axis="x").points = [1, 2, 4, 10]
        msg = "Currently unable to ensure isotropic filtering."
        with self.assertRaisesRegex(RuntimeError, msg):
            raymond(self.source, epsilon=1, isotropic=True)

    def test_allclose_call(self):
        # Test that the tolerance from the config is used in the isotropy checks.
        stub_config = copy.copy(ants.config.GlobalConfiguration())
        stub_config.__init__()
        stub_config["ants_tolerance"]["raymond_filter_isotropy_tolerance"] = 0.1
        config_patch = mock.patch("ants.config.CONFIG", new=stub_config)
        with config_patch:
            with mock.patch("ants.utils.ndarray.allclose") as mock_allclose:
                raymond(self.source, epsilon=1, isotropic=True)
        mock_allclose.assert_called_with(60.0, [60.0], tolerance=0.1)


if __name__ == "__main__":
    ants.tests.main()
