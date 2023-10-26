# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import datetime
import os
import tempfile
import unittest.mock as mock

import ants.io.save as save
import ants.tests
import iris.tests.stock as stock


class TestAll(ants.tests.TestCase):
    def setUp(self):
        self.filename = "my_filename"
        patch = mock.patch("iris.save")
        self.netcdf = patch.start()
        self.addCleanup(patch.stop)

        self.cubes = [ants.tests.stock.geodetic((3, 2))]

    def test_default_args_with_no_nc_extension(self):
        save.netcdf(self.cubes, self.filename)
        self.netcdf.assert_called_once_with(
            self.cubes,
            f"{self.filename}.nc",
            saver="nc",
            netcdf_format="NETCDF4_CLASSIC",
            local_keys=None,
            unlimited_dimensions=None,
            zlib=False,
            complevel=4,
            fill_value=None,
        )

    def test_default_args_with_nc_extension(self):
        filename = f"{self.filename}.nc"
        save.netcdf(self.cubes, filename)
        self.netcdf.assert_called_once_with(
            self.cubes,
            filename,
            saver="nc",
            netcdf_format="NETCDF4_CLASSIC",
            local_keys=None,
            unlimited_dimensions=None,
            zlib=False,
            complevel=4,
            fill_value=None,
        )

    def test_update_history_False_with_no_history_message(self):
        save.netcdf(self.cubes, self.filename, update_history=False)
        self.assertNotIn("history", self.cubes[0].attributes)

    def test_update_history_False_with_history_message(self):
        history_message = "My history message"
        with mock.patch("warnings.warn") as mock_warn:
            save.netcdf(
                self.cubes,
                self.filename,
                update_history=False,
                history_message=history_message,
            )
        self.assertEqual(mock_warn.call_count, 1)
        self.assertNotIn("history", self.cubes[0].attributes)

    def test_update_history_True_with_no_history_message(self):
        arguments = "program arg1 arg2"
        with mock.patch("sys.argv", new=arguments.split()):
            # 'update_history' is 'True' by default.
            save.netcdf(self.cubes, self.filename)
        output = self.cubes[0].attributes["history"]
        # Ignore the timestamp at the start of the 'history' attribute:
        msg = f"End of {output} != {arguments}"
        self.assertTrue(output.endswith(arguments), msg)

    def test_update_history_True_with_history_message(self):
        history_message = "My history message"
        # 'update_history' is 'True' by default.
        save.netcdf(self.cubes, self.filename, history_message=history_message)
        output = self.cubes[0].attributes["history"]
        # Ignore the timestamp at the start of the 'history' attribute:
        msg = f"End of {output} != {history_message}"
        self.assertTrue(output.endswith(history_message), msg)

    def test_netcdf_classic_coerce_call(self):
        # Ensure that a cube being written with NETCDF4_CLASSIC calls the
        # coerce function
        patch_func = "ants.io.save._coerce_netcdf_classic_dtypes"
        with mock.patch(patch_func) as patched:
            save.netcdf(self.cubes, self.filename)
        self.assertTrue(patched.called)

    def test_netcdf_classic_workaround_call(self):
        # Ensure that a cube being written with NETCDF4_CLASSIC calls the
        # iris workaround function
        patch_func = "ants.io.save._iris_netcdf4_classic_workaround"
        with mock.patch(patch_func) as patched:
            save.netcdf(self.cubes, self.filename)
        self.assertTrue(patched.called)


class TestSave(ants.tests.TestCase):
    def setUp(self):
        self.cube = stock.lat_lon_cube()
        self.cdl_filename_to_delete = None

    def tearDown(self):
        # If a test is performing a CDL comparison of the netCDF created as a
        # side effect of the ancillary saving, we cannot use
        # NamedTemporaryFile to automatically delete the netCDF file.  In this
        # case, self.cdl_filename_to_delete should be set to the netCDF
        # filename to ensures that the temporary file is deleted after the
        # test, regardless of whether the test passes or fails.
        if self.cdl_filename_to_delete:
            os.remove(self.cdl_filename_to_delete)

    def test_save_netcdf(self):
        fh = tempfile.NamedTemporaryFile(suffix=".nc")

        # Fix the date so that we can test
        today = datetime.datetime(1, 1, 1)
        with mock.patch("ants.fileformats.datetime") as dt_patch:
            dt_patch.today.return_value = today
            with mock.patch("sys.argv", new=["program", "arg1", "arg2"]):
                save.netcdf(self.cube, fh.name)
        self.assertCDL(fh.name, ("fileformats", "save_netcdf.cdl"))

    def test_fill_value_set_correctly_in_netcdf_file(self):
        fh = tempfile.NamedTemporaryFile(suffix=".nc")

        # Fix the date so that we can test
        today = datetime.datetime(1, 1, 1)
        with mock.patch("ants.fileformats.datetime") as dt_patch:
            dt_patch.today.return_value = today
            with mock.patch("sys.argv", new=["program", "arg1", "arg2"]):
                self.cube.attributes["STASH"] = "m01s01i001"
                save.netcdf(self.cube, fh.name, fill_value=-99)
        self.assertCDL(fh.name, ("fileformats", "fill_value_set_netcdf.cdl"))


if __name__ == "__main__":
    ants.tests.main()
