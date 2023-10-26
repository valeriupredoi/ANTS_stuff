# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
from ants.analysis._raymond import filter_length_scale_to_eps_iso


class TestAll(ants.tests.TestCase):
    def test_value(self):
        filter_length_scale = 3
        delta_lambda = 0.83333333
        earth_radius = 6371229.0
        lats = [89.58333333, 47.91666667, 6.25, -35.41666667, -77.08333333]
        eps = filter_length_scale_to_eps_iso(
            filter_length_scale, lats, delta_lambda, earth_radius
        )
        tar = [1.68075e-10, 1.67772160e05, 8.51735512e02, 4.83460911e-06, 1.37062240e03]
        self.assertArrayAlmostEqual(eps, tar, decimal=5)


if __name__ == "__main__":
    ants.tests.main()
