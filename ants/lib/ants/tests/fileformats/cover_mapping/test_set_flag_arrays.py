# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import unittest.mock as mock

import ants.tests
import numpy as np
from ants.fileformats.cover_mapping import set_flag_arrays


class TestFlagValues(ants.tests.TestCase):
    def setUp(self):
        self.cube = mock.Mock(name="cube", attributes=dict())
        self.target = [1, 2, 3, 4]

    def assert_flag_arrays_call(self, flag_values):
        flag_meanings = "one two three four"
        set_flag_arrays(self.cube, flag_values, flag_meanings)
        self.assertEqual(self.cube.attributes["flag_values"], self.target)

    def test_list_type(self):
        flag_values = [1, 2, 3, 4]
        self.assert_flag_arrays_call(flag_values)

    def test_numpy_type(self):
        flag_values = np.array([1, 2, 3, 4])
        self.assert_flag_arrays_call(flag_values)

    def test_float_type(self):
        flag_values = np.array([1.0, 2.0, 3.0, 4.0])
        self.assert_flag_arrays_call(flag_values)

    def test_float_not_whole_number(self):
        flag_values = np.array([1.5, 2.6, 3.9, 4.1])
        msg = "Flag values cannot be represented as integers"
        with self.assertRaisesRegex(ValueError, msg):
            self.assert_flag_arrays_call(flag_values)


class TestFlagMeanings(ants.tests.TestCase):
    def setUp(self):
        self.cube = mock.Mock(name="cube", attributes=dict())
        self.target = "one two three four"

    def assert_flag_arrays_call(self, flag_meanings):
        flag_values = [1, 2, 3, 4]
        set_flag_arrays(self.cube, flag_values, flag_meanings)
        self.assertEqual(self.cube.attributes["flag_meanings"], self.target)

    def test_string_type(self):
        flag_values = "one two three four"
        self.assert_flag_arrays_call(flag_values)

    def test_redundant_spaces(self):
        flag_values = "     one             two   three      four"
        self.assert_flag_arrays_call(flag_values)

    def test_iterable_type(self):
        flag_values = ["one", "two", "three", "four"]
        self.assert_flag_arrays_call(flag_values)

    def test_alt_delimiter(self):
        flag_values = "one,two,three,four"
        self.assert_flag_arrays_call(flag_values)

    def test_alt_delimiter_with_space(self):
        flag_values = "one, two, three, four"
        self.assert_flag_arrays_call(flag_values)

    def test_more_than_one_delimiter(self):
        flag_values = "one,two three;four"

        msg = "More than one potential delimiter found"
        with self.assertRaisesRegex(ValueError, msg):
            self.assert_flag_arrays_call(flag_values)


class TestExceptions(ants.tests.TestCase):
    def setUp(self):
        self.cube = mock.Mock(name="cube", attributes=dict())

    def test_meanings_values_differ_length(self):
        flag_values = [1, 2, 3, 4, 5]
        flag_meanings = "one two three four"

        msg = "Flag values and flag meanings are not the same length"
        with self.assertRaisesRegex(ValueError, msg):
            set_flag_arrays(self.cube, flag_values, flag_meanings)


if __name__ == "__main__":
    ants.tests.main()
