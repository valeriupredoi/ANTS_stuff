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

        # Move around in a circle
        for direction in range(1, 9):
            location.shift(direction)

        # The final location is where we started. At one point we moved off
        # the edge of the domain but as j_overlap is used in the LocaionClass
        # then the program manages to retain its original location.
        target_i = 4
        target_j = 4

        self.assertEqual(location.i, target_i)
        self.assertEqual(location.j, target_j)


if __name__ == "__main__":
    ants.tests.main()
