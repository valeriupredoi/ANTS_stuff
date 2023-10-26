# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
import iris
import numpy as np
from proc_ants.order_nemo_rivers import get_min_distance


class TestAll(ants.tests.TestCase):
    def test_all(self):

        # Generate a river going through middle of array
        river_sequence = np.array(
            [
                [0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0],
                [1, 2, 3, 4, 5],
                [0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0],
            ]
        )

        sequence_cube = iris.cube.Cube(
            river_sequence, long_name="river sequence number"
        )
        scaling = 1.0
        j = 2
        i = 4

        min_distance = get_min_distance(j, i, sequence_cube, scaling)

        self.assertAlmostEqual(min_distance, 8)


if __name__ == "__main__":
    ants.tests.main()


if __name__ == "__main__":
    ants.tests.main()
