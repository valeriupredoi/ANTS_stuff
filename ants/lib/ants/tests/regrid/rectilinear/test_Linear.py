# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import unittest.mock as mock

import ants.tests
from ants.regrid.rectilinear import Linear


class Test___init__(ants.tests.TestCase):
    def test_default_args(self):
        scheme = Linear()
        self.assertEqual(scheme.extrapolation_mode, "nan")


class Test_regridder(ants.tests.TestCase):
    def test_call_args(self):
        src = mock.sentinel.src_cube
        tgt = mock.sentinel.tgt_cube
        patch_ants_regridder = mock.patch(
            "ants.regrid.rectilinear._RectilinearRegridder"
        )

        scheme = Linear()
        with patch_ants_regridder as mock_ants_regridder:
            scheme.regridder(src, tgt)
        mock_ants_regridder.assert_called_once_with(src, tgt, "linear", "nan")


class Test_interpolator(ants.tests.TestCase):
    def test_call_args(self):
        src = mock.Mock("Cube")
        src.coord = mock.Mock(side_effect=lambda x: x)
        coords = [mock.sentinel.coord1, mock.sentinel.coord2]
        patch_ants_interpolator = mock.patch(
            "ants.regrid.rectilinear._RectilinearInterpolator"
        )

        scheme = Linear()
        with patch_ants_interpolator as mock_ants_interpolator:
            scheme.interpolator(src, coords)

        mock_ants_interpolator.assert_called_once_with(src, coords, "linear", "nan")


class Test___repr__(ants.tests.TestCase):
    def test(self):
        scheme = Linear()
        tar = "Linear(extrapolation_mode=nan)"
        self.assertEqual(repr(scheme), tar)


if __name__ == "__main__":
    ants.tests.main()
