# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import unittest.mock as mock
from copy import copy

import ants.tests
import numpy as np
from ants.analysis.cover_mapping import CoverMapper


class _Common(object):
    def setUp(self):
        self.source_types = ["desert", "vegetation", "sea"]
        self.target_types = [1, 101]
        self.cover_map = np.array([[100, 0], [100, 0], [0, 100]])
        self.cover_mapper = CoverMapper(
            self.source_types, self.target_types, self.cover_map
        )


class Test___init__(_Common, ants.tests.TestCase):
    def test_missing_source(self):
        # Ensure that we fall over if there are less source types than there
        # are mappings.
        cover_map = self.cover_map[1:]
        msg = (
            "The cover map has 2 source mappings to 2 target "
            "classifications while there are 3 expected source types and 2 "
            "expected target types."
        )
        with self.assertRaisesRegex(RuntimeError, msg):
            CoverMapper(self.source_types, self.target_types, cover_map)

    def test_missing_target(self):
        # Ensure that we fall over if there are less target types than there
        # are mappings.
        cover_map = self.cover_map[:, 1:]
        msg = (
            "The cover map has 3 source mappings to 1 target "
            "classifications while there are 3 expected source types and 2 "
            "expected target types."
        )
        with self.assertRaisesRegex(RuntimeError, msg):
            CoverMapper(self.source_types, self.target_types, cover_map)

    def test_invalid_cover_map(self):
        # Ensure that we fall over if the mappings do not add to a 100%
        # contribution for each target type.
        # Here we mix fractional with percentage contributions.
        cover_map = self.cover_map.astype("float")
        cover_map[0, 0] = 0.5
        msg = (
            "Incorrectly defined cover map table, contributions do not "
            "add up to 100%"
        )
        with self.assertRaisesRegex(RuntimeError, msg):
            CoverMapper(self.source_types, self.target_types, cover_map)

    def test_percentages_cover_map(self):
        # Assumed percentages.
        cm = CoverMapper(self.source_types, self.target_types, self.cover_map)
        self.assertArrayEqual(cm.cover_map, self.cover_map)

    def test_fractionals_cover_map(self):
        # Assumed fractionals.
        cover_map = self.cover_map.copy()
        cover_map = cover_map / 100.0
        cm = CoverMapper(self.source_types, self.target_types, cover_map)
        self.assertArrayEqual(cm.cover_map, self.cover_map)

    def test_fractionals_rounding(self):
        # Ensure that converted percentages add up to 100.
        cover_map = self.cover_map.astype("float")
        cover_map = cover_map / 100.0
        cover_map[0, 0] = 0.981
        cover_map[0, 1] = 0.019
        cm = CoverMapper(self.source_types, self.target_types, cover_map)
        target_cover_map = [[98, 2], [100, 0], [0, 100]]
        self.assertArrayEqual(cm.cover_map, target_cover_map)

    def test_error_if_non_int_targets(self):
        # Ensure an exception is thrown if the target types are not integer
        target_types = self.target_types
        target_types[0] = "sea"
        with self.assertRaises(ValueError):
            CoverMapper(self.source_types, target_types, self.cover_map)


class Test_source_types(_Common, ants.tests.TestCase):
    def test_all(self):
        self.assertArrayEqual(self.cover_mapper.source_types, self.source_types)


class Test_target_types(_Common, ants.tests.TestCase):
    def test_all(self):
        self.assertArrayEqual(self.cover_mapper.target_types, self.target_types)


class Test_cover_map(_Common, ants.tests.TestCase):
    def test_all(self):
        self.assertArrayEqual(self.cover_mapper.cover_map, self.cover_map)


class Test_reorder_cover_map(_Common, ants.tests.TestCase):
    def test_reordering(self):
        types_target = np.roll(self.source_types, 1)
        covermap_target = np.roll(self.cover_map, 1, axis=0)
        self.cover_mapper.reorder_cover_map(types_target)
        self.assertArrayEqual(self.cover_mapper.source_types, types_target)
        self.assertArrayEqual(self.cover_mapper.cover_map, covermap_target)

    def test_nochange(self):
        types_target = copy(self.source_types)
        covermap_target = copy(self.cover_map)
        self.cover_mapper.reorder_cover_map(types_target)
        self.assertArrayEqual(self.cover_mapper.source_types, types_target)
        self.assertArrayEqual(self.cover_mapper.cover_map, covermap_target)

    def test_underspecified_source_types(self):
        # Not all source types were specified.
        msg = "Underspecified source types"
        with self.assertRaisesRegex(ValueError, msg):
            self.cover_mapper.reorder_cover_map(self.source_types[0:2])

    def test_unexpected_source_types(self):
        # Not all the provided source types match those in the existing
        # CoverMap object.
        source_types = copy(self.source_types)
        source_types[0] = "moon"
        msg = "Source types do not match cover table description:"
        with self.assertRaisesRegex(ValueError, msg):
            self.cover_mapper.reorder_cover_map(source_types)


class Test_apply_transform(ants.tests.TestCase):
    def setUp(self):
        self.data = np.array([[1, 1, 3], [9, 9, 3], [1, 3, 1]], dtype="int8")
        self.flag_values = np.array([1, 3, 9], dtype="int8")

    @staticmethod
    def apply_transform(data, cover_mapping, flag_values):
        source_types_stub = [str(i) for i in range(cover_mapping.shape[0])]
        target_types_stub = list(range(cover_mapping.shape[1]))
        cover_mapper = CoverMapper(source_types_stub, target_types_stub, cover_mapping)
        return cover_mapper.apply_transform(data, flag_values)

    def assert_apply_transform_dtype(
        self, flag_values_dtype=None, crosswalk_dtype=None
    ):
        # Apply the transform for the case where we have 1-to-1 mapping.
        cover_map = np.array([[100, 0], [100, 0], [0, 100]], dtype="int8")
        if flag_values_dtype is not None:
            self.flag_values = self.flag_values.astype(flag_values_dtype)
        if crosswalk_dtype is not None:
            cover_map = cover_map.astype(crosswalk_dtype)

        res = self.apply_transform(self.data, cover_map, self.flag_values)
        tgt = np.array(
            [
                [[100, 100, 100], [0, 0, 100], [100, 100, 100]],
                [[0, 0, 0], [100, 100, 0], [0, 0, 0]],
            ]
        )
        self.assertArrayEqual(res, tgt)

    def test_value(self):
        self.assert_apply_transform_dtype()

    def test_dtype_crosswalk(self):
        self.assert_apply_transform_dtype(crosswalk_dtype="int64")

    def test_dtype_flag_values_dtype(self):
        self.assert_apply_transform_dtype(flag_values_dtype="int64")

    def test_value_alt_flag_value_order(self):
        flag_values = np.array([9, 1, 3], dtype="int8")
        cover_map = np.array([[0, 100], [100, 0], [100, 0]], dtype="int8")
        res = self.apply_transform(self.data, cover_map, flag_values)
        tgt = np.array(
            [
                [[100, 100, 100], [0, 0, 100], [100, 100, 100]],
                [[0, 0, 0], [100, 100, 0], [0, 0, 0]],
            ]
        )
        self.assertArrayEqual(res, tgt)

    def test_single_src_multiple_target_mapping(self):
        cover_map = np.array([[100, 0], [50, 50], [0, 100]], dtype="int8")
        res = self.apply_transform(self.data, cover_map, self.flag_values)
        tgt = np.array(
            [
                [[100, 100, 50], [0, 0, 50], [100, 50, 100]],
                [[0, 0, 50], [100, 100, 50], [0, 50, 0]],
            ]
        )
        self.assertArrayEqual(res, tgt)

    def test_masked_array_retention(self):
        # Ensure that masked elements remain on the array returned.
        cover_map = np.array([[100, 0], [100, 0], [0, 100]], dtype="int8")
        self.data = np.ma.array(self.data)
        self.data[:] = np.ma.masked

        res = self.apply_transform(self.data, cover_map, self.flag_values)
        tgt = np.ma.array(
            [
                [[100, 100, 100], [0, 0, 100], [100, 100, 100]],
                [[0, 0, 0], [100, 100, 0], [0, 0, 0]],
            ]
        )
        tgt[:] = np.ma.masked
        self.assertMaskedArrayEqual(res, tgt)


class Test_save(_Common, ants.tests.TestCase):
    def testall(self):
        jpatch = mock.patch("json.dump")
        opatch = mock.patch("ants.analysis.cover_mapping.open", mock.mock_open())

        path = mock.sentinel.transform_path
        with jpatch as json_patch, opatch as open_patch:
            self.cover_mapper.save(path)
        open_patch.assert_called_once_with(path, "w")
        args = json_patch.call_args_list[0][0]

        tar = (
            {
                "source": ["desert", "vegetation", "sea"],
                "target": [1, 101],
                "cover_map": [[100, 0], [100, 0], [0, 100]],
            },
            open_patch().__enter__(),
        )
        self.assertEqual(args, tar)


class Test_load(_Common, ants.tests.TestCase):
    def test_valid(self):
        transform = {
            "source": self.source_types,
            "target": self.target_types,
            "cover_map": self.cover_map,
        }
        jpatch = mock.patch("json.load", return_value=transform)
        opatch = mock.patch("ants.analysis.cover_mapping.open", mock.mock_open())

        path = mock.sentinel.transform_path
        with jpatch, opatch as open_patch:
            res = self.cover_mapper.load(path)

        open_patch.assert_called_once_with(path, "r")
        self.assertArrayEqual(res.source_types, self.source_types)
        self.assertArrayEqual(res.target_types, self.target_types)
        self.assertArrayEqual(res.cover_map, self.cover_map)

    def test_invalid(self):
        transform = {}
        jpatch = mock.patch("json.load", return_value=transform)
        opatch = mock.patch("ants.analysis.cover_mapping.open", mock.mock_open())

        path = mock.sentinel.transform_path
        with jpatch, opatch:
            msg = "source. Valid covermap files should include"
            with self.assertRaisesRegex(KeyError, msg):
                self.cover_mapper.load(path)


if __name__ == "__main__":
    ants.tests.main()
