# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import os
import tempfile

import ants.io.save as save
import ants.tests
import ants.tests.stock as stock
import check_cumf
import iris


class TestAll(ants.tests.TestCase):
    def setUp(self):
        self.cube1 = stock.geodetic((2, 3))
        self.cube1.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi(
            "m01s00i274"
        )
        self.cube2 = self.cube1.copy()

        self.tmpdir1 = tempfile.TemporaryDirectory()
        self.tmpdir2 = tempfile.TemporaryDirectory()
        self.filename = tempfile.NamedTemporaryFile(dir=self.tmpdir1.name)
        self.filename.close()
        self.ref_filepath = self.tmpdir2.name

    def check_cumf(self):
        save.ancil(self.cube1, self.filename.name)
        filename2 = os.path.join(
            self.ref_filepath, os.path.basename(self.filename.name)
        )
        save.ancil(self.cube2, filename2)
        check_cumf.main(self.filename.name, filename2)

    def test_nodiff(self):
        self.check_cumf()

    def test_diff(self):
        self.cube1.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi("m01s00i27")
        msg = "FAIL: Ancillaries do not match:"
        with self.assertRaisesRegex(RuntimeError, msg):
            self.check_cumf()


if __name__ == "__main__":
    ants.tests.main()
