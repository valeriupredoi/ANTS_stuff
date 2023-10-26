# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import os
import shutil
import tempfile
import unittest.mock as mock

import ants.fileformats.ancil
import ants.io.save
import ants.tests
from ants.tests.io.old_save.common import ancil_save_call, run_command


class TestNewSaveAncil(ants.tests.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.cubes = [ants.tests.stock.geodetic((3, 2), stash="m12s34i567")]

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def _save_using_old_saver(self, filename):
        old_filename = os.path.join(self.test_dir, filename)
        with mock.patch("warnings.warn") as mock_warn:
            ants.fileformats.ancil.save(self.cubes, old_filename)
        self.assertEqual(mock_warn.call_count, 1)
        calls = [ancil_save_call()]
        mock_warn.assert_has_calls(calls)
        return old_filename

    def _save_using_new_saver(self, filename):
        new_filename = os.path.join(self.test_dir, filename)
        ants.io.save.ancil(self.cubes, new_filename)
        return new_filename

    def test_save_ancil_with_extension(self):
        old_filename = self._save_using_old_saver("old_ancil_with_ext.ancil")
        new_filename = self._save_using_new_saver("new_ancil_with_ext.ancil")
        compare_ancil(old_filename, new_filename)

    def test_save_ancil_without_extension(self):
        old_filename = self._save_using_old_saver("old_ancil_without_ext")
        new_filename = self._save_using_new_saver("new_ancil_without_ext")
        compare_ancil(old_filename, new_filename)


def compare_ancil(file1, file2):
    command = ["mule-cumf", file1, file2]
    run_command(command, AssertionError)


if __name__ == "__main__":
    ants.tests.main()
