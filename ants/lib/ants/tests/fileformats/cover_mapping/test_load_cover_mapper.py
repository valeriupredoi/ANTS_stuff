# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import unittest.mock as mock

import ants.tests
from ants.fileformats.cover_mapping import load_cover_mapper


class TestAll(ants.tests.TestCase):
    def setUp(self):
        self.transform = {
            "source": ["desert", "vegetation", "sea"],
            "target": [1, 101],
            "cover_map": [[0, 100], [100, 0], [100, 0]],
        }

        patch = mock.patch("json.load", return_value=self.transform)
        self.mock_json = patch.start()
        self.addCleanup(patch.stop)

        patch = mock.patch("builtins.open")
        self.mock_open = patch.start()
        self.addCleanup(patch.stop)

    def test_all(self):
        path = mock.sentinel.transform_path
        with mock.patch("ants.analysis.cover_mapping.CoverMapper") as pcm:
            load_cover_mapper(path)
        pcm.load.assert_called_once_with(path)


if __name__ == "__main__":
    ants.tests.main()
