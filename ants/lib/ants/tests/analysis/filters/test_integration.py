# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
import numpy as np
from ants.analysis.filters import raymond


class TestRaymond(ants.tests.TestCase):
    def setUp(self):
        self.epsilon = 1.0
        length = 50
        shift = 0
        self.xvals = np.arange(shift, length + shift)
        self.perturbations = self.wave(25, length)
        self.bare = self.wave(3, length)
        self.bare = self.bare.reshape(1, -1)
        self.data = self.bare + self.perturbations
        self.source = ants.tests.stock.geodetic(data=self.data)

    @staticmethod
    def wave(number, length):
        # Note: endpoint False as bc's are written with u_0 = u_(n+1), not
        # u_0 = u_n
        return 100 * np.cos(np.linspace(0, number * 2 * np.pi, length, endpoint=False))

    def test_values(self):
        # Ensure that when periodic is True that x-filter is periodic.
        raymond(self.source, 1)
        self.assertTrue((abs(self.data - self.bare) < 1e-2).all())


if __name__ == "__main__":
    ants.tests.main()
