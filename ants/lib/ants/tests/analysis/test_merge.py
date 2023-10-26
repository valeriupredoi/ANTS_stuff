# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import unittest.mock as mock

import ants.tests
import iris
from ants.analysis import merge


class Testall(ants.tests.TestCase):
    def setUp(self):
        patch = mock.patch(
            "ants.analysis._merge.merge", return_value=mock.sentinel.merged_cube
        )
        self.mock_merge = patch.start()
        self.addCleanup(patch.stop)

        patch = mock.patch("ants.analysis.FillMissingPoints")
        self.mock_fill = patch.start()
        self.addCleanup(patch.stop)

        patch = mock.patch(
            "ants.utils.cube.sort_cubes", side_effect=lambda x, y: (x, y)
        )
        self.mock_sort = patch.start()
        self.addCleanup(patch.stop)

    def test_call_args(self):
        primary_cube = mock.Mock(name="primary_cube", spec=iris.cube.Cube)
        alternate_cube = mock.Mock(name="alternate_cube", spec=iris.cube.Cube)
        res = merge(primary_cube, alternate_cube, mock.sentinel.validity_polygon)

        self.assertFalse(self.mock_fill.called)
        self.assertIs(res, mock.sentinel.merged_cube)
        self.mock_merge.assert_called_once_with(
            primary_cube, alternate_cube, mock.sentinel.validity_polygon
        )


if __name__ == "__main__":
    ants.tests.main()
