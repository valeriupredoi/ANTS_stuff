# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
"""
Stub doc for testing.
"""
import argparse
import unittest.mock as mock

import ants.tests
from ants.command_parse import AntsArgParser
from ants.tests.io.old_save.common import saver_call


class TestAntsArgParser(ants.tests.TestCase):
    def setUp(self):
        patch = mock.patch("ants.config.CONFIG.parse_configuration")
        self.mock_config = patch.start()
        self.addCleanup(patch.stop)

        patch = mock.patch("ants.config.filepath_readable")
        self.mock_filepath_readable = patch.start()
        self.mock_filepath_readable.side_effect = lambda x: x
        self.addCleanup(patch.stop)

        patch = mock.patch("ants.config.dirpath_writeable")
        self.mock_dirpath_writeable = patch.start()
        self.mock_dirpath_writeable.side_effect = lambda x: x
        self.addCleanup(patch.stop)

    def test_saver_deprecation(self):
        new = [
            "program",
            "/path/to/source",
            "--target-lsm",
            "/path/to/lsm",
            "-o",
            "/path/to/output",
            "--saver",
            "ukca",
        ]

        with mock.patch("sys.argv", new=new):
            parser = AntsArgParser(target_lsm=True)
            with mock.patch("warnings.warn") as mock_warn:
                args = parser.parse_args()
            self.assertEqual(mock_warn.call_count, 1)
            calls = [saver_call()]
            mock_warn.assert_has_calls(calls)

        target_args = argparse.Namespace(
            ants_config=None,
            land_threshold=None,
            target_lsm="/path/to/lsm",
            output="/path/to/output",
            saver="ukca",
            sources=["/path/to/source"],
            use_new_saver=False,
            netcdf_only=False,
        )
        self.assertFalse(self.mock_config.called)
        self.assertEqual(args, target_args)


if __name__ == "__main__":
    ants.tests.main()
