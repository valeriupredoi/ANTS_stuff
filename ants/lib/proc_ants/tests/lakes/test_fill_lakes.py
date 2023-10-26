# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
from proc_ants.lakes import fill_lakes


class TestAll(ants.tests.TestCase):
    def setUp(self):
        self.cube = ants.tests.stock.geodetic((10, 10), xlim=(28, 38), ylim=(-6, 5))
        self.cube.data = self.cube.data.astype("int8")
        self.cube.data[:] = 2
        self.cube.attributes["flag_meanings"] = "ocean lake other"
        self.cube.attributes["flag_values"] = [1, 2, 3]

    def test_fill_lake_all_lake(self):
        # Ensure that the entire geometry is filled when all values are
        # lake - i.e. all set to ocean.
        seed = [-1, 33]
        res = fill_lakes(self.cube, "victoria", seed, "lake", "ocean", constrain=True)
        tar = [
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
            [2, 2, 2, 2, 1, 1, 2, 2, 2, 2],
            [2, 2, 2, 2, 1, 1, 2, 2, 2, 2],
            [2, 2, 2, 2, 1, 1, 2, 2, 2, 2],
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
        ]
        self.assertArrayEqual(res.data, tar)

    def test_fill_lake_partial_lake(self):
        # Ensure that only lake values are changed, leaving other values
        # within the geometry unchanged.
        self.cube.data[1::2] = 3
        seed = [-1, 33]
        res = fill_lakes(self.cube, "victoria", seed, "lake", "ocean", constrain=True)
        tar = [
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
            [2, 2, 2, 2, 1, 1, 2, 2, 2, 2],
            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
        ]
        self.assertArrayEqual(res.data, tar)

    def test_inapropriate_seed(self):
        seed = [-90, 0]
        msg = r".*is not contained within the"
        with self.assertRaisesRegex(ValueError, msg):
            fill_lakes(self.cube, "victoria", seed, "lake", "ocean")


if __name__ == "__main__":
    ants.tests.main()
