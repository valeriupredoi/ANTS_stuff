# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
import numpy as np
from proc_ants.location_class import LocationClass


class TestAll(ants.tests.TestCase):
    def test_all(self):

        orca_mask_array = np.array(
            [
                [0, 0, 0, 0, 0],
                [1, 0, 0, 0, 0],
                [1, 1, 0, 0, 0],
                [1, 1, 1, 0, 0],
                [1, 1, 1, 1, 0],
            ]
        )

        start_i = 2
        start_j = 2
        location = LocationClass(start_j, start_i, orca_mask_array.shape)

        sea_coord_list = location.find_neighbouring_sea(orca_mask_array, 2)
        sea_coord_array = np.array(sea_coord_list)

        target_sea_coord_array = np.array(
            [[3, 2], [2, 1], [3, 1], [4, 3], [4, 2], [1, 0], [2, 0]]
        )

        self.assertArrayEqual(sea_coord_array, target_sea_coord_array)


if __name__ == "__main__":
    ants.tests.main()
