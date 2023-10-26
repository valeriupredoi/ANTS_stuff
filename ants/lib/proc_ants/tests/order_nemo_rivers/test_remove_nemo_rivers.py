# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
import iris
import numpy as np
from proc_ants.order_nemo_rivers import remove_nemo_rivers


class TestAll(ants.tests.TestCase):
    def test_all(self):

        river_index_array = np.array(
            [
                [4, 0, 0, 0, 0],
                [4, 2, 2, 0, 0],
                [4, 2, 2, 0, 0],
                [4, 5, 6, 0, 0],
                [4, 5, 6, 0, 0],
            ]
        )

        river_number = np.array(
            [
                [4, 0, 0, 0, 0],
                [4, 2, 2, 0, 0],
                [4, 2, 2, 0, 0],
                [4, 0, 6, 0, 0],
                [4, 0, 6, 0, 0],
            ]
        )

        river_number_cube = iris.cube.Cube(river_number, long_name="river number")

        remove_nemo_rivers(river_index_array, river_number_cube)

        target_river_number_array = np.array(
            [
                [4, 0, 0, 0, 0],
                [4, 2, 2, 0, 0],
                [4, 2, 2, 0, 0],
                [4, 0, 5, 0, 0],
                [4, 0, 5, 0, 0],
            ]
        )

        self.assertArrayEqual(river_index_array, target_river_number_array)
        self.assertArrayEqual(river_number_cube.data, target_river_number_array)


if __name__ == "__main__":
    ants.tests.main()


if __name__ == "__main__":
    ants.tests.main()
