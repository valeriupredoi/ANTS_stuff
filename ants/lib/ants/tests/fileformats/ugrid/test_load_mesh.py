# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
from ants.fileformats._ugrid import load_mesh


class TestAll(ants.tests.TestCase):
    def test_exception_for_multiple_filenames(self):
        with self.assertRaisesRegex(ValueError, "Multiple mesh"):
            load_mesh(["foo", "bar"])


if __name__ == "__main__":
    ants.tests.main()
