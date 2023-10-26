# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
import numpy as np
from proc_ants.order_nemo_rivers import assign_river_number_to_points


class TestAll(ants.tests.TestCase):
    def test_all(self):

        river_index_array = np.zeros([7, 7])

        runoff_nc_array = np.array(
            [
                [4, 0, 0, 0, 0, 0, 0],
                [4, 2, 2, 0, 0, 0, 0],
                [4, 2, 2, 0, 0, 0, 1],
                [4, 5, 6, 0, 0, 1, 1],
                [4, 5, 6, 0, 0, 1, 1],
                [4, 5, 3, 3, 3, 1, 1],
                [4, 3, 3, 3, 1, 1, 1],
            ]
        )

        twod_index = (5, 3)

        river_number = 9

        amount = 10

        accumulative_amount_list = []

        assign_river_number_to_points(
            river_index_array,
            runoff_nc_array,
            twod_index,
            river_number,
            amount,
            accumulative_amount_list,
        )

        target_river_index_array = np.array(
            [
                [9, 0, 0, 0, 0, 0, 0],
                [9, 9, 9, 0, 0, 0, 0],
                [9, 9, 9, 0, 0, 0, 9],
                [9, 9, 9, 0, 0, 9, 9],
                [9, 9, 9, 0, 0, 9, 9],
                [9, 9, 9, 9, 9, 9, 9],
                [9, 9, 9, 9, 9, 9, 9],
            ]
        )

        self.assertArrayEqual(river_index_array, target_river_index_array)
        self.assertArrayEqual(accumulative_amount_list[0], [9, 98])


if __name__ == "__main__":
    ants.tests.main()
