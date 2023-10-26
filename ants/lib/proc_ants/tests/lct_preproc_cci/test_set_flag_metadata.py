# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import unittest.mock as mock

import ants.tests
from proc_ants.lct_preproc_cci import set_flag_metadata


class TestAll(ants.tests.TestCase):
    def testall(self):
        cube = mock.Mock("cube", attributes={})
        set_flag_metadata(cube)
        self.assertIn("flag_meanings", cube.attributes)
        self.assertIn("flag_values", cube.attributes)


if __name__ == "__main__":
    ants.tests.main()
