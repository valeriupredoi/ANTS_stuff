# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.

import ants.tests
import ants.tests.stock as stock
from ants.analysis import stdev


class Test_stdev(ants.tests.TestCase):
    def test_negative_variance(self):
        # Issue particularly presents itself due to the squaring of the 32bit
        # array then subsequently being regridded into a 64 bit array.
        source = stock.geodetic((5, 5))
        source.data.astype("float32")
        mean = stock.geodetic((5, 5))
        res = stdev(source, mean)
        self.assertFalse((res.data < 0).any())


if __name__ == "__main__":
    ants.tests.main()
