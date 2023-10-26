# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
import numpy as np
from proc_ants.location_class import LocationClass


class TestAll(ants.tests.TestCase):
    def test_all(self):

        river_index_array = np.array(
            [
                [1, 0, 0, 0, 0, 0, 0],
                [4, 2, 2, 0, 0, 0, 0],
                [4, 2, 2, 0, 0, 0, 0],
                [4, 5, 6, 0, 0, 0, 0],
                [4, 5, 6, 0, 0, 0, 0],
                [4, 5, 0, 0, 0, 0, 0],
                [4, 0, 0, 0, 0, 0, 0],
            ]
        )

        # 1 = sea, 0 = land
        orca_mask_array = np.array(
            [
                [1, 0, 0, 0, 0, 0, 0],
                [1, 1, 1, 0, 0, 0, 0],
                [1, 1, 1, 0, 0, 0, 0],
                [1, 1, 1, 0, 0, 0, 0],
                [1, 1, 1, 0, 0, 0, 0],
                [1, 1, 1, 1, 1, 0, 0],
                [1, 1, 1, 1, 1, 0, 0],
            ]
        )

        start_i = 4
        start_j = 6
        location = LocationClass(start_j, start_i, orca_mask_array.shape)

        river_number = location.diamond_search(
            river_index_array=river_index_array,
            orca_mask_array=orca_mask_array,
            min_distance=2,
        )

        # The target river number is 2 as the search finds the minimum
        # river number except that river number 1 is too far away
        # The diamond search looks two grid boxes away plus the find_neighbouring_sea
        # function looks a further 3 grid boxes away (total = 5 grid boxes)
        target_river_number = 2

        self.assertEqual(river_number, target_river_number)


if __name__ == "__main__":
    ants.tests.main()
