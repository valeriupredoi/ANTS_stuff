# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import unittest.mock as mock

import ants.tests
from ants.fileformats.cover_mapping import get_flag_arrays


class TestValue(ants.tests.TestCase):
    def testall(self):
        attributes = {
            "flag_meanings": "one two Three four Five",
            "flag_values": [1, 2, 3, 4, 5],
        }
        cube = mock.Mock(name="cube", attributes=attributes)
        flg_val, flg_mng = get_flag_arrays(cube)
        self.assertArrayEqual(flg_mng, ["one", "two", "three", "four", "five"])
        self.assertArrayEqual(flg_val, [1, 2, 3, 4, 5])


class TestExceptions(ants.tests.TestCase):
    def test_meanings_values_differ_length(self):
        attributes = {
            "flag_meanings": "one two three four",
            "flag_values": [1, 2, 3, 4, 5],
        }
        cube = mock.Mock(name="cube", attributes=attributes)

        msg = "Missing flag value/meaning pair"
        with self.assertRaisesRegex(RuntimeError, msg):
            get_flag_arrays(cube)


if __name__ == "__main__":
    ants.tests.main()
