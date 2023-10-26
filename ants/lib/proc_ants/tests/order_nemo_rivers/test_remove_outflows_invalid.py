# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
import numpy as np
from proc_ants.order_nemo_rivers import remove_outflows_invalid


class TestAll(ants.tests.TestCase):
    def test_all(self):

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

        remove_outflows_invalid(None, runoff_nc_array)

        target_runoff_nc_array = np.array(
            [
                [0, 0, 0, 0, 0, 0, 0],
                [0, 2, 2, 0, 0, 0, 0],
                [0, 2, 2, 0, 0, 0, 0],
                [0, 5, 6, 0, 0, 1, 0],
                [0, 5, 6, 0, 0, 1, 0],
                [0, 5, 3, 3, 3, 1, 0],
                [0, 0, 0, 0, 0, 0, 0],
            ]
        )

        self.assertArrayEqual(runoff_nc_array, target_runoff_nc_array)


if __name__ == "__main__":
    ants.tests.main()
