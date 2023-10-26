# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
import numpy as np
from ants.analysis._raymond import epsilon_iso


class TestAll(ants.tests.TestCase):
    def test_value(self):
        # We pick this shape semi arbitrarily except for the fact that it
        # produces value of epsilon which are suitable for testing (i.e. not
        # max eps for every latitude).
        shape = (21, 43)
        epsilon = 1.0
        eps = epsilon_iso(shape, epsilon)

        # We check values of eps as it goes towards the poles.
        # Values, were originally derived by comparison with the CAP, using a
        # 216x432 array (see #759).
        tar = np.array(
            [
                1.67772160e05,
                3.81124295e04,
                1.75824693e03,
                2.31520287e02,
                5.14933201e01,
                1.59596391e01,
                6.34246185e00,
                3.11919924e00,
                1.87631622e00,
                1.38009795e00,
                1.24511718e00,
                1.38009795e00,
                1.87631622e00,
                3.11919924e00,
                6.34246185e00,
                1.59596391e01,
                5.14933201e01,
                2.31520287e02,
                1.75824693e03,
                3.81124295e04,
                1.67772160e05,
            ]
        )
        self.assertArrayAlmostEqual(eps, tar, decimal=4)


if __name__ == "__main__":
    ants.tests.main()
