# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import unittest.mock as mock

import ants.tests
import iris
import numpy as np
from ants.analysis._merge import _UGridFillMissingPoints


class Common(object):
    def setUp(self):
        self.source = ants.tests.stock.mesh_C4(
            load_data=True, constraint="sample_data"
        )[0]
        source_mask = np.array([False] * 96)
        self.source.data = np.ma.array(self.source.data, mask=source_mask)


class TestNoTargetMask(Common, ants.tests.TestCase):
    def test_no_mask_and_nan_cell(self):
        """need to test nans converted to mask"""
        expected = ants.tests.stock.mesh_C4(load_data=True, constraint="sample_data")[0]
        expected.data[25] = 22.0

        self.source.data = self.source.data.data
        self.source.data[25] = np.nan

        filler = _UGridFillMissingPoints(self.source)
        filler(self.source)
        actual = self.source

        self.assertMaskedArrayEqual(expected.data, actual.data)

    def test_no_mask_and_no_nan_cells(self):
        self.source.data = self.source.data.data
        expected = hash(self.source)

        filler = _UGridFillMissingPoints(self.source)
        filler(self.source)
        actual = hash(self.source)

        self.assertEqual(expected, actual)

    def test_single_masked_cell(self):
        expected = ants.tests.stock.mesh_C4(load_data=True, constraint="sample_data")[0]
        expected.data[10] = 7.0

        self.source.data.mask[10] = True

        filler = _UGridFillMissingPoints(self.source)
        filler(self.source)
        actual = self.source

        self.assertMaskedArrayEqual(expected.data, actual.data)

    def test_same_cell_masked_and_nan(self):
        expected = ants.tests.stock.mesh_C4(load_data=True, constraint="sample_data")[0]
        # Will be defined as NaN:
        expected.data[0] = 65.0
        # Will be masked and NaN
        expected.data[95] = 48.0

        self.source.data[0] = np.nan
        self.source.data[95] = np.nan
        self.source.data.mask[95] = True
        filler = _UGridFillMissingPoints(self.source)
        filler(self.source)
        actual = self.source

        self.assertMaskedArrayEqual(expected.data, actual.data)

    def test_query_different_cube(self):
        expected = ants.tests.stock.mesh_C4(load_data=True, constraint="sample_data")[0]
        expected.data = expected.data + 100
        expected.data[10] = 107

        self.source.data.mask[10] = True
        filler = _UGridFillMissingPoints(self.source)
        new_source = ants.tests.stock.mesh_C4(load_data=True, constraint="sample_data")[
            0
        ]
        new_source.data = np.ma.array(new_source.data + 100, mask=self.source.data.mask)
        filler(new_source)
        actual = new_source

        self.assertMaskedArrayEqual(expected.data, actual.data)

    def test_ndsupport(self):
        expected = np.ma.arange(1, 97)
        expected[30] = 30
        expected = np.tile(expected, 2).reshape(2, 96)
        expected[1, :] += 20

        self.source.data[30] = np.ma.masked
        source2 = self.source.copy()
        self.source.add_aux_coord(iris.coords.AuxCoord(0, long_name="bing"), None)
        source2.add_aux_coord(iris.coords.AuxCoord(2, long_name="bing"), None)
        source2.data += 20
        source = iris.cube.CubeList([self.source, source2]).merge_cube()
        filler = _UGridFillMissingPoints(source)
        filler(source)
        actual = source.data
        self.assertMaskedArrayEqual(expected, actual)


class TestTargetMaskWith1DSource(Common, ants.tests.TestCase):
    def test_target_and_source_masked(self):
        """
        Cells masked in both source and target should remain masked.

        """
        expected = ants.tests.stock.mesh_C4(load_data=True, constraint="sample_data")[0]
        expected.data = np.ma.masked_array(expected.data)

        target_mask = expected.copy(data=np.zeros_like(expected.data, dtype="bool"))
        target_mask.data[20] = True
        self.source.data[20] = np.ma.masked
        expected.data[20] = np.ma.masked

        filler = _UGridFillMissingPoints(self.source, target_mask)
        filler(self.source)
        actual = self.source

        self.assertMaskedArrayEqual(expected.data, actual.data)

    def test_target_unmasked_and_source_masked(self):
        """
        Cells masked in source and unmasked in target should be filled and
        unmasked.

        """
        expected = ants.tests.stock.mesh_C4(load_data=True, constraint="sample_data")[0]
        expected.data = np.ma.masked_array(expected.data)

        target_mask = expected.copy(data=np.zeros_like(expected.data, dtype="bool"))
        self.source.data[20] = np.ma.masked
        expected.data[20] = 25

        filler = _UGridFillMissingPoints(self.source, target_mask)
        filler(self.source)
        actual = self.source

        self.assertMaskedArrayEqual(expected.data, actual.data)

    def test_target_masked_and_source_unmasked(self):
        """
        Cells unmasked in source and masked in target should be masked.

        """
        expected = ants.tests.stock.mesh_C4(load_data=True, constraint="sample_data")[0]
        expected.data = np.ma.masked_array(expected.data)

        target_mask = expected.copy(data=np.zeros_like(expected.data, dtype="bool"))
        target_mask.data[20] = True
        expected.data[20] = np.ma.masked

        filler = _UGridFillMissingPoints(self.source, target_mask)
        filler(self.source)
        actual = self.source

        self.assertMaskedArrayEqual(expected.data, actual.data)

    def test_target_and_source_unmasked(self):
        """
        Cells unmasked in source and unmasked should remain unmasked.

        """
        expected = ants.tests.stock.mesh_C4(load_data=True, constraint="sample_data")[0]
        expected.data = np.ma.masked_array(expected.data)

        target_mask = expected.copy(data=np.zeros_like(expected.data, dtype="bool"))

        filler = _UGridFillMissingPoints(self.source, target_mask)
        filler(self.source)
        actual = self.source

        self.assertMaskedArrayEqual(expected.data, actual.data)


class TestTargetMaskWithNDSource(ants.tests.TestCase):
    """
    Target mask should be 1D and applied to every step of the extra
    dimensions.

    """

    def setUp(self):
        source = ants.tests.stock.mesh_C4()
        source.data = np.arange(1, 97)
        self.target_mask = ants.tests.stock.mesh_C4()
        self.target_mask.data = np.zeros_like(source.data, dtype="bool")
        source2 = source.copy()
        source.add_aux_coord(iris.coords.AuxCoord(0, long_name="bing"), None)
        source2.add_aux_coord(iris.coords.AuxCoord(2, long_name="bing"), None)
        source2.data += 20
        self.source = iris.cube.CubeList([source, source2]).merge_cube()
        self.source.data = np.ma.array(self.source.data)

    def test_target_and_source_masked(self):
        """
        Cells masked in both source and target should remain masked.

        """
        expected = self.source.copy()
        self.target_mask.data[20] = True
        self.source.data[:, 20] = np.ma.masked
        expected.data[:, 20] = np.ma.masked

        filler = _UGridFillMissingPoints(self.source, self.target_mask)
        filler(self.source)
        actual = self.source

        self.assertMaskedArrayEqual(expected.data, actual.data)

    def test_target_unmasked_and_source_masked(self):
        """
        Cells masked in source and unmasked in target should be filled and
        unmasked.

        """
        expected = self.source.copy()
        expected.data = np.ma.masked_array(expected.data)

        self.source.data[:, 20] = np.ma.masked
        expected.data[0, 20] = 25
        expected.data[1, 20] = 45

        filler = _UGridFillMissingPoints(self.source, self.target_mask)
        filler(self.source)
        actual = self.source

        self.assertMaskedArrayEqual(expected.data, actual.data)

    def test_target_masked_and_source_unmasked(self):
        """
        Cells unmasked in source and masked in target should be masked.

        """
        expected = self.source.copy()
        expected.data = np.ma.masked_array(expected.data)

        self.target_mask.data[20] = True
        expected.data[:, 20] = np.ma.masked

        filler = _UGridFillMissingPoints(self.source, self.target_mask)
        filler(self.source)
        actual = self.source

        self.assertMaskedArrayEqual(expected.data, actual.data)

    def test_target_and_source_unmasked(self):
        """
        Cells unmasked in source and unmasked should remain unmasked.

        """
        expected = self.source.copy()
        expected.data = np.ma.masked_array(expected.data)

        filler = _UGridFillMissingPoints(self.source, self.target_mask)
        filler(self.source)
        actual = self.source

        self.assertMaskedArrayEqual(expected.data, actual.data)


class TestInit(Common, ants.tests.TestCase):
    @mock.patch("ants.analysis._merge._UGridFillMissingPoints._create_kdtree")
    def test_create_tree_called_by_init(self, patch_kdtree):
        _UGridFillMissingPoints(self.source)
        patch_kdtree.assert_called_once()

    def test_accepts_target_mask(self):
        _UGridFillMissingPoints(
            self.source,
            target_mask=self.source.copy(data=np.zeros_like(self.source.data)),
        )

    @mock.patch("ants.analysis._merge._UGridFillMissingPoints._create_kdtree")
    @mock.patch("ants.utils.cube.is_equal_hgrid", return_value=False)
    def test_rejects_different_horizontal_spec_target_mask(self, *args):
        target_mask = np.ma.masked_array(np.zeros_like(self.source))
        with self.assertRaisesRegex(ValueError, "The provided target_mask a"):
            _UGridFillMissingPoints(self.source, target_mask=target_mask)

    # Wanting to test how target mask is used - so mock out return from cartopy
    # since the lat/lon to Cartesian transform adds unnecessary complexity to
    # the test otherwise:
    @mock.patch(
        "ants.analysis._merge._UGridFillMissingPoints." "_transform_points",
        return_value=(np.arange(1, 97).repeat(3).reshape(96, 3)),
    )
    def test_target_mask_not_added_to_masked_points(self, _):
        expected = np.ma.arange(1, 97).repeat(3).reshape(96, 3)

        # Target mask contains a single cell that is no masked in source or
        # expected.  Source has a single masked cell that is masked in
        # expected:
        target_mask = self.source.copy(np.zeros_like(self.source.data))
        target_mask.data[2] = 1
        self.source.data[10] = np.ma.masked
        expected[10, :] = np.ma.masked
        expected = expected[~expected.mask].reshape(-1, 3).data

        with mock.patch("ants.analysis._merge.KDTree") as patch_KDTree:
            _UGridFillMissingPoints(self.source, target_mask)
        # Can't use assert_called_once_with due to array ambiguities - hence
        # split into 3 asserts for equivalent functionality:
        patch_KDTree.assert_called_once()
        self.assertArrayEqual(patch_KDTree.call_args[0][0], expected)
        self.assertEqual(1, len(patch_KDTree.call_args[0]))


class TestCreateKDTree(Common, ants.tests.TestCase):
    def test_kdtree_not_called_when_not_required(self):
        self.source.data = self.source.data.data
        with mock.patch("ants.analysis._merge.KDTree") as kdpatch:
            _UGridFillMissingPoints(self.source)
        self.assertFalse(kdpatch.called)

    @mock.patch("ants.analysis._merge.KDTree")
    @mock.patch(
        "ants.analysis._merge._UGridFillMissingPoints." "_transform_points",
        return_value=np.ma.masked_array(np.arange(96 * 3).reshape((96, 3))),
    )
    def test_KDTree_call_single_dimension(self, _, patch_KDTree):
        expected = np.arange(96 * 3).reshape(96, 3)
        expected = np.ma.array(np.vstack((expected[:10], expected[11:])))

        self.source.data[10] = np.ma.masked
        _UGridFillMissingPoints(self.source)

        # Can't use assert_called_once_with due to array ambiguities - hence
        # split into 3 asserts for equivalent functionality:
        patch_KDTree.assert_called_once()
        self.assertArrayEqual(patch_KDTree.call_args[0][0], expected)
        self.assertEqual(1, len(patch_KDTree.call_args[0]))

    @mock.patch("ants.analysis._merge.KDTree")
    @mock.patch(
        "ants.analysis._merge._UGridFillMissingPoints." "_transform_points",
        return_value=np.ma.masked_array(np.arange(96 * 3).reshape((96, 3))),
    )
    def test_KDTree_call_multidimensional(self, _, patch_KDTree):
        expected = np.arange(96 * 3).reshape(96, 3)
        expected = np.ma.array(np.vstack((expected[:10], expected[11:])))

        self.source.data[10] = np.ma.masked
        source2 = self.source.copy()
        self.source.add_aux_coord(iris.coords.AuxCoord(0, long_name="bing"), None)
        source2.add_aux_coord(iris.coords.AuxCoord(2, long_name="bing"), None)
        source = iris.cube.CubeList([self.source, source2]).merge_cube()
        _UGridFillMissingPoints(source)

        # Can't use assert_called_once_with due to array ambiguities - hence
        # split into 3 asserts for equivalent functionality:
        patch_KDTree.assert_called_once()
        self.assertArrayEqual(patch_KDTree.call_args[0][0], expected)
        self.assertEqual(1, len(patch_KDTree.call_args[0]))


class TestExceptions(Common, ants.tests.TestCase):
    def test_query_different_cube_coordinates(self):
        filler = _UGridFillMissingPoints(self.source)
        new_source = ants.tests.stock.mesh_C4(load_data=True, constraint="sample_data")[
            0
        ]
        points = new_source.coord("longitude").points + 0.1
        new_source.coord("longitude").points = points
        with self.assertRaisesRegex(ValueError, "Coordinates differ "):
            filler(new_source)

    def test_not_ugrid_data(self):
        source = ants.tests.stock.geodetic((2, 2))
        with self.assertRaisesRegex(ValueError, "Source not recognised"):
            _UGridFillMissingPoints(source)

    def test_no_valid_data(self):
        self.source.data[:] = np.ma.masked
        with self.assertRaisesRegex(ValueError, "No valid data provided for"):
            _UGridFillMissingPoints(self.source)

    def test_inconsistent_masks(self):
        new_source = ants.tests.stock.mesh_C4(load_data=True, constraint="sample_data")[
            0
        ]
        new_source_mask = np.array([False] * 96)
        new_source_mask[25] = True
        new_source.data = np.ma.array(new_source.data, mask=new_source_mask)

        filler = _UGridFillMissingPoints(self.source)

        with self.assertRaisesRegex(ValueError, "Masks differ "):
            filler(new_source)


if __name__ == "__main__":
    ants.tests.main()
