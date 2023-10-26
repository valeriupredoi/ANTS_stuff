# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import unittest.mock as mock

import ants.tests
import iris
import numpy as np
from ants.analysis.cover_mapping import CoverMapper, SCTTransformer


class Test___init__(ants.tests.TestCase):
    def testall(self):
        transform = mock.sentinel.transform
        sct_transformer = SCTTransformer(transform)
        self.assertIs(sct_transformer._transform, transform)


class Test___call__(ants.tests.TestCase):
    def setUp(self):
        # Source cube
        self.src_cube = ants.tests.stock.geodetic((3, 3))
        self.src_cube.data = np.array([[1, 1, 3], [9, 9, 3], [1, 3, 1]], dtype="int8")
        # CF specified that flag_meanings be space separated
        self.src_cube.attributes["flag_meanings"] = "desert vegetation sea"
        self.src_cube.attributes["flag_values"] = [1, 3, 9]

        # Source transformation
        source_types = np.array(["desert", "vegetation", "sea"])
        target_types = np.array([5, 10])
        cover_map = np.array([[100, 0], [100, 0], [0, 100]])
        self.transform = CoverMapper(source_types, target_types, cover_map)

        # Target grid cube
        self.target_cube = mock.Mock(name="target_grid_cube")
        # Classifications cube - from splitting the source cube into their
        # relevant types.
        self.classifications = iris.cube.Cube(
            np.array(
                [[[1, 1, 1], [0, 0, 1], [1, 1, 1]], [[0, 0, 0], [1, 1, 0], [0, 0, 0]]]
            )
        )
        class_coord = iris.coords.AuxCoord(
            ["desert", "sea"], long_name="pseudo_level", units=1
        )
        self.classifications.add_aux_coord(class_coord, 0)
        self.classifications.add_dim_coord(self.src_cube.coord(axis="x"), 2)
        self.classifications.add_dim_coord(self.src_cube.coord(axis="y"), 1)

        patch = mock.patch("ants.analysis.mean")
        self.mock_mean = patch.start()
        self.mock_mean.side_effect = lambda x, y: x
        self.addCleanup(patch.stop)

        self.sct_transformer = SCTTransformer(self.transform)

    def assertCoverMap(self, res):
        target_data = [
            [[1, 1, 1], [0, 0, 1], [1, 1, 1]],
            [[0, 0, 0], [1, 1, 0], [0, 0, 0]],
        ]
        target_classes = [5, 10]
        self.assertArrayEqual(res.data, target_data)
        self.assertArrayEqual(res.coord("pseudo_level").points, target_classes)

    def test_cube_metadata(self):
        res = self.sct_transformer(self.src_cube, self.target_cube)
        self.assertCoverMap(res)

    def test_mean_call(self):
        self.sct_transformer(self.src_cube, self.target_cube)
        self.assertTrue(self.mock_mean.called)

    def _mean_cube(self, data):
        data = np.asarray(data)
        cube = ants.tests.stock.geodetic(data=data)
        pseudo_level_coord = iris.coords.DimCoord(
            range(data.shape[0]), long_name="pseudo_level"
        )
        cube.add_dim_coord(pseudo_level_coord, 0)
        return cube

    def test_normalisation(self):
        # Ensure that values are normalised by adding to non-zero
        # classifications and in the case where there are none leave alone.
        mean_return = self._mean_cube(
            [
                [[0, 90, 0], [100, 45, 100], [0, 0, 0]],
                [[100, 0, 90], [0, 45, 0], [0, 0, 0]],
            ]
        )
        target = [
            [[0, 1, 0], [1, 0.5, 1], [0, 0, 0]],
            [[1, 0, 1], [0, 0.5, 0], [0, 0, 0]],
        ]
        with mock.patch("ants.analysis.mean", return_value=mean_return):
            with mock.patch("warnings.warn") as mock_warn:
                res = self.sct_transformer(self.src_cube, self.target_cube)
        self.assertArrayAlmostEqual(res.data, target)
        msg = (
            "Locations present with no classification fraction, ignoring "
            "such locations."
        )
        mock_warn.assert_called_once_with(msg)

    def test_dtype_boolean_special_case(self):
        # Ensure that where the crosswalk allows and no target grid is defined
        # that the result is of int8 type.  Real examples include the
        # pre-processing of the landseamask.
        res = self.sct_transformer(self.src_cube)
        self.assertEqual(res.data.dtype, np.dtype("int8"))

    def test_mask_removal(self):
        self.src_cube.data = np.ma.array(self.src_cube.data, mask=False)
        res = self.sct_transformer(self.src_cube)
        self.assertFalse(np.ma.isMaskedArray(res.data))


if __name__ == "__main__":
    ants.tests.main()
