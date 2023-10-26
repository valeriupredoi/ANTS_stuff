# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
import numpy as np
from ants.analysis.cover_mapping import CoverMapper, SCTTransformer


class TestAll(ants.tests.TestCase):
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
        cover_map = np.array([[1, 0], [1, 0], [0, 1]])
        self.transform = CoverMapper(source_types, target_types, cover_map)

    def _call_class_transform(self):
        tar = [[[1, 1, 1], [0, 0, 1], [1, 1, 1]], [[0, 0, 0], [1, 1, 0], [0, 0, 0]]]
        sct_transformer = SCTTransformer(self.transform)
        res = sct_transformer(self.src_cube)
        self.assertArrayEqual(res.data, tar)

    def test_grouped_mappings(self):
        self._call_class_transform()

    def test_reordered_mappings(self):
        # Test the case where the cube has different flag ordering to the
        # crosswalk table.
        self.src_cube.attributes["flag_meanings"] = "sea desert vegetation"
        self.src_cube.attributes["flag_values"] = [9, 1, 3]
        self._call_class_transform()


if __name__ == "__main__":
    ants.tests.main()
