# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
from ants.analysis._raymond import filter_length_scale_to_eps


class TestAll(ants.tests.TestCase):
    def test_value(self):
        filter_length_scale = 1
        delta_lambda = 0.83333333
        earth_radius = 6371229.0
        eps = filter_length_scale_to_eps(
            filter_length_scale, delta_lambda, earth_radius
        )
        tar = 2.0635708098788648
        self.assertEqual(eps, tar)


if __name__ == "__main__":
    ants.tests.main()
