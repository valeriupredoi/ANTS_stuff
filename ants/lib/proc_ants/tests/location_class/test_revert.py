# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
from proc_ants.location_class import LocationClass


class TestAll(ants.tests.TestCase):
    def test_all(self):

        # Set an initial location near the top of a grid
        location = LocationClass(4, 4, [6, 6])

        # Go north 4 grid boxes, going off the edge of the domain
        location.north(amount=4)

        # Go east
        location.east()

        # Check we have moved
        target_i = 5
        target_j = 5
        target_j_overlap = 8

        self.assertEqual(location.i, target_i)
        self.assertEqual(location.j, target_j)
        self.assertEqual(location.j_overlap, target_j_overlap)

        location.revert()

        # The final location is where we started as we reverted
        target_i = 4
        target_j = 4
        target_j_overlap = 4

        self.assertEqual(location.i, target_i)
        self.assertEqual(location.j, target_j)
        self.assertEqual(location.j_overlap, target_j_overlap)


if __name__ == "__main__":
    ants.tests.main()
