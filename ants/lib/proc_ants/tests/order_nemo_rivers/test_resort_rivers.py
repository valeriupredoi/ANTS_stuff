# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
import numpy as np
from proc_ants.order_nemo_rivers import resort_rivers


class TestAll(ants.tests.TestCase):
    def test_all(self):

        river_index_array = np.array(
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

        accumulative_amount_list = [
            [3, 0.5],
            [1, 0.2],
            [5, 0.3],
            [6, 0.1],
            [2, 0.15],
            [4, 0.9],
        ]

        river_index_new_array = resort_rivers(
            river_index_array, accumulative_amount_list
        )

        target_river_index_array = np.array(
            [
                [1, 0, 0, 0, 0, 0, 0],
                [1, 5, 5, 0, 0, 0, 0],
                [1, 5, 5, 0, 0, 0, 4],
                [1, 3, 6, 0, 0, 4, 4],
                [1, 3, 6, 0, 0, 4, 4],
                [1, 3, 2, 2, 2, 4, 4],
                [1, 2, 2, 2, 4, 4, 4],
            ]
        )

        self.assertArrayEqual(river_index_new_array, target_river_index_array)


if __name__ == "__main__":
    ants.tests.main()
