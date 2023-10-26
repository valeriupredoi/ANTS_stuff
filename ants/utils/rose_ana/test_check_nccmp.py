# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import tempfile

import ants.io.save as save
import ants.tests
import ants.tests.stock as stock
import check_nccmp


class TestAll(ants.tests.TestCase):
    def setUp(self):
        self.cube1 = stock.geodetic((2, 3))
        self.cube2 = self.cube1.copy()
        self.tmpdir1 = tempfile.TemporaryDirectory()
        self.filename1 = tempfile.NamedTemporaryFile(
            dir=self.tmpdir1.name, suffix=".nc"
        )
        self.filename1.close()
        self.tmpdir2 = tempfile.TemporaryDirectory()
        self.filename2 = tempfile.NamedTemporaryFile(
            dir=self.tmpdir2.name, suffix=".nc"
        )
        self.filename2.close()

    def test_nodiff(self):
        save.netcdf(self.cube1, self.filename1.name)
        save.netcdf(self.cube2, self.filename2.name)
        assert check_nccmp.main(self.filename1.name, self.filename2.name, None) is None

    def test_diff(self):
        self.cube1.rename("dummy_name")
        save.netcdf(self.cube1, self.filename1.name)
        save.netcdf(self.cube2, self.filename2.name)
        msg = "NetCDF comparison failure."
        with self.assertRaisesRegex(RuntimeError, msg):
            check_nccmp.main(self.filename1.name, self.filename2.name, None)

    def test_specified_attributes_raises_error_when_not_ignored(self):
        self.cube1.attributes["fake_attribute"] = "dummy"
        save.netcdf(self.cube1, self.filename1.name)
        save.netcdf(self.cube2, self.filename2.name)
        with self.assertRaisesRegex(RuntimeError, "NetCDF comparison failure."):
            check_nccmp.main(self.filename1.name, self.filename2.name, None)

    def test_ignore_specified_attributes(self):
        self.cube1.attributes["fake_attribute"] = "dummy"
        attributes = ["fake_attribute"]
        save.netcdf(self.cube1, self.filename1.name)
        save.netcdf(self.cube2, self.filename2.name)
        assert (
            check_nccmp.main(self.filename1.name, self.filename2.name, attributes)
            is None
        )

    def test_history_attribute_ignored_by_nccmp_by_default(self):
        self.cube1.attributes["history"] = "dummy"
        save.netcdf(self.cube1, self.filename1.name)
        save.netcdf(self.cube2, self.filename2.name)
        assert check_nccmp.main(self.filename1.name, self.filename2.name, None) is None


if __name__ == "__main__":
    ants.tests.main()
