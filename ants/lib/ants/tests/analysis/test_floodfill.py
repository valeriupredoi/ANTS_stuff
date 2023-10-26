# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
import numpy as np
from ants.analysis import floodfill
from ants.exceptions import FloodfillError


class TestValues(ants.tests.TestCase):
    def testall(self):
        # Ensure that the array that is returned has the floodfill applied to
        # it.
        arr = np.ones((5, 4))
        arr[1:3, 1:3] = 10

        tar = arr.copy()
        tar[tar == 10] = 5

        floodfill(arr, (2, 2), 5)
        self.assertArrayEqual(arr, tar)

    def test_nofill(self):
        # Raise a suitable exception where the starting location does not
        # appear to require filling (i.e. already has the fillvalue).
        arr = np.ones((5, 4))
        arr[2, 2] = 5
        msg = "The value at location"
        with self.assertRaisesRegex(FloodfillError, msg):
            floodfill(arr, (2, 2), 5)


if __name__ == "__main__":
    ants.tests.main()
