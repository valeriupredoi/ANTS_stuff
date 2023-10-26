# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
import numpy as np
from proc_ants.order_nemo_rivers import scale_by_resolution


class TestAll(ants.tests.TestCase):
    def test_all(self):
        # Array of latitudes
        nav_lat = np.array(
            [
                [0.2, 0.2, 0.2, 0.2, 0.2],
                [0.1, 0.1, 0.1, 0.1, 0.1],
                [0.0, 0.0, 0.0, 0.0, 0.0],
                [-0.1, -0.1, -0.1, -0.1, -0.1],
                [-0.2, -0.2, -0.2, -0.2, -0.2],
            ]
        )

        nav_lon = np.array(
            [
                [89.8, 89.9, 90.0, 90.1, 90.2],
                [89.8, 89.9, 90.0, 90.1, 90.2],
                [89.8, 89.9, 90.0, 90.1, 90.2],
                [89.8, 89.9, 90.0, 90.1, 90.2],
                [89.8, 89.9, 90.0, 90.1, 90.2],
            ]
        )

        scaling = scale_by_resolution(nav_lat, nav_lon)

        self.assertAlmostEqual(scaling, 2.5)


if __name__ == "__main__":
    ants.tests.main()
