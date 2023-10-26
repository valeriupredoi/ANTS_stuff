# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
import numpy as np
from ants.analysis._raymond import filters_non_periodic as filters


class Test_filters(ants.tests.TestCase):
    # Checking the return values of filters where the target values were taken
    # from the CAP in the case of non periodic case.
    def setUp(self):
        self.arr = np.array(
            [144, 147, 164, 167, 167, 109, 183, 121, 136, 187], dtype="float"
        )

    def test_values(self):
        res = filters(self.arr, 1)
        tar = np.array([0, 14, -39, 416, -1014, 1642, -1674, 254, 36, 0])
        self.assertArrayEqual(res, tar)


if __name__ == "__main__":
    ants.tests.main()
