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

        is_coastal_array = np.empty_like(orca_mask_array)

        for j in range(orca_mask_array.shape[0]):
            for i in range(orca_mask_array.shape[1]):
                location = LocationClass(j, i, orca_mask_array.shape)
                is_coastal_array[j, i] = location.is_coastal(orca_mask_array)

        target_is_coastal_array = np.array(
            [
                [1, 1, 0, 0, 1],
                [0, 1, 1, 0, 1],
                [0, 0, 1, 1, 1],
                [0, 0, 0, 1, 1],
                [0, 0, 0, 0, 1],
            ]
        )

        self.assertArrayEqual(is_coastal_array, target_is_coastal_array)


if __name__ == "__main__":
    ants.tests.main()
