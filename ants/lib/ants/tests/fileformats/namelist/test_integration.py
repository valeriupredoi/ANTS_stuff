# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import unittest.mock as mock
from io import StringIO

import ants.tests
from ants.fileformats.namelist import (
    load_cap_horizontal,
    load_lfric_vertical,
    load_um_vertical,
)


@ants.tests.skip_f90nml
class TestHorizontal(ants.tests.TestCase):
    def run_regular(self, callback=None):
        data = """&GRID
 POINTS_LAMBDA_TARG=2,POINTS_PHI_TARG=2,PHI_ORIGIN_TARG=45
/
"""
        file1 = StringIO(data)
        filename = "dummy_filename"
        patch_open = mock.patch("f90nml.parser.open", create=True, return_value=file1)
        with patch_open:
            cube = next(load_cap_horizontal(filename, callback=callback))
        return cube

    def test_regular(self):
        cube = self.run_regular()
        self.assertCML(
            cube, ("fileformats", "namelist", "integration_regular.cml"), checksum=False
        )

    def test_callback(self):
        def my_callback(cube, groups, filename):
            cube.attributes["dummy_change"] = "flop"

        cube = self.run_regular(callback=my_callback)
        self.assertEqual(cube.attributes["dummy_change"], "flop")

    def test_variable_res(self):
        data = """&GRID
 POINTS_LAMBDA_TARG=2,POINTS_PHI_TARG=2
/
"""
        file1 = StringIO(data)

        data = """&HORIZGRID
 LAMBDA_INPUT_P=1, 2, 3
 LAMBDA_INPUT_U=1.5, 2.5, 3.5
 PHI_INPUT_P=1, 2, 3
 PHI_INPUT_V=0.5, 1.5, 2.5
/
"""
        file2 = StringIO(data)

        filenames = ["filename1", "filename2"]
        patch_open = mock.patch(
            "f90nml.parser.open", create=True, side_effect=[file1, file2]
        )
        with patch_open:
            cube = next(load_cap_horizontal(filenames))

        self.assertCML(
            cube,
            ("fileformats", "namelist", "integration_variable_res.cml"),
            checksum=False,
        )


@ants.tests.skip_f90nml
class TestVertical(ants.tests.TestCase):
    def setUp(self, callback=None):
        data = """&VERTLEVS
z_top_of_model = 5.0
first_constant_r_rho_level = 3
eta_theta = 0.0, 0.1, 0.225, 0.4, 0.6, 1.0
eta_rho = 0.05, 0.1625, 0.3125, 0.5, 0.8
/
"""
        self.file1 = StringIO(data)
        self.filename = "filename"

    def _run_um_vertical(self, callback=None):
        patch_open = mock.patch(
            "f90nml.parser.open", create=True, return_value=self.file1
        )
        with patch_open:
            cube = next(load_um_vertical(self.filename, callback=callback))
        return cube

    def _run_lfric_vertical(self, callback=None):
        patch_open = mock.patch(
            "f90nml.parser.open", create=True, return_value=self.file1
        )
        with patch_open:
            cube = next(load_lfric_vertical(self.filename, callback=callback))
        return cube

    def test_um_vertical(self):
        cube = self._run_um_vertical()
        self.assertCML(
            cube,
            ("fileformats", "namelist", "integration_um_vertical.cml"),
            checksum=False,
        )

    def test_lfric_vertical(self):
        cube = self._run_lfric_vertical()
        self.assertCML(
            cube,
            ("fileformats", "namelist", "integration_lfric_vertical.cml"),
            checksum=False,
        )

    def test_callback(self):
        def my_callback(cube, groups, filename):
            cube.attributes["dummy_change"] = "flop"

        cube = self._run_um_vertical(callback=my_callback)
        self.assertEqual(cube.attributes["dummy_change"], "flop")


if __name__ == "__main__":
    ants.tests.main()
