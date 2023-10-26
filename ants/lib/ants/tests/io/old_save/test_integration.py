# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import datetime
import os
import tempfile
import unittest.mock as mock

import ants
import ants.tests
import ants.tests.stock
import iris
import iris.tests.stock as stock
from ants.tests.io.old_save.common import (
    ancil_save_call,
    ants_save_call,
    netcdf_save_call,
    ukca_netcdf_save_call,
)


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

    def test_cf_netcdf_save(self):
        fh = tempfile.NamedTemporaryFile(suffix=".nc")

        # Fix the date so that we can test
        today = datetime.datetime(1, 1, 1)
        with mock.patch("ants.fileformats.datetime") as dt_patch:
            dt_patch.today.return_value = today
            with mock.patch("sys.argv", new=["program", "arg1", "arg2"]):
                with mock.patch("warnings.warn") as mock_warn:
                    ants.save(self.cube, fh.name)
                self.assertEqual(mock_warn.call_count, 2)
                calls = [ants_save_call(), netcdf_save_call()]
                mock_warn.assert_has_calls(calls)
        self.assertCDL(fh.name, ("fileformats", "save_netcdf.cdl"))

    def test_cf_netcdf_save_with_broken_ancil(self):
        with mock.patch("ants.fileformats.ancil.save", side_effect=IOError):
            with mock.patch("ants.fileformats.netcdf.cf.save") as nc_patch:
                # Ignore the tautology - it's to trap the error to allow test
                # to proceed to the real assertion.
                with self.assertRaises(IOError):
                    cube = self.cube
                    cube.attributes["STASH"] = "m01s01i001"
                    with mock.patch("warnings.warn") as mock_warn:
                        ants.save(cube, "dummy_fnme", saver="ancil")
                    # 'ancil.save' and 'netcdf.save' are mocked:
                    self.assertEqual(mock_warn.call_count, 1)
                    calls = [ants_save_call()]
                    mock_warn.assert_has_calls(calls)
                nc_patch.assert_called_once_with(self.cube, "dummy_fnme.nc")

    def test_ancil_save_ancil_args(self):
        # Ensure expected arguments passed to the underlying um.save
        apatch = mock.patch("ants.fileformats.ancil.save")
        cfpatch = mock.patch("ants.fileformats.netcdf.cf.save")
        with apatch as ancil_patch, cfpatch:
            filename = "dummy_fnme"
            cube = self.cube
            cube.attributes["STASH"] = "m01s01i001"
            with mock.patch("warnings.warn") as mock_warn:
                ants.save(cube, filename, saver="ancil")
            # 'ancil.save' and 'netcdf.save' are mocked:
            self.assertEqual(mock_warn.call_count, 1)
            calls = [ants_save_call()]
            mock_warn.assert_has_calls(calls)
            ancil_patch.assert_called_once_with(self.cube, filename)

    def test_ukca_netcdf_save(self):
        ukcapatch = mock.patch("ants.fileformats.netcdf.ukca._ukca_conventions")
        cfpatch = mock.patch("ants.fileformats.netcdf.cf.save")
        with cfpatch, ukcapatch as nc_patch:
            with mock.patch("warnings.warn") as mock_warn:
                ants.save(self.cube, "dummy_fnme", saver="ukca")
            # 'netcdf.save' is mocked:
            self.assertEqual(mock_warn.call_count, 2)
            calls = [ants_save_call(), ukca_netcdf_save_call()]
            mock_warn.assert_has_calls(calls)
        nc_patch.assert_called_once()

    def test_no_saver_specified(self):
        cfpatch = mock.patch("ants.fileformats.netcdf.cf.save")
        apatch = mock.patch("ants.fileformats.ancil.save")
        with apatch as ancil_patch, cfpatch:
            cube = self.cube
            cube.attributes["STASH"] = "m01s01i001"
            with mock.patch("warnings.warn") as mock_warn:
                ants.save(cube, "dummy_fnme")
            self.assertEqual(mock_warn.call_count, 1)
            calls = [ants_save_call()]
            mock_warn.assert_has_calls(calls)
        ancil_patch.assert_called_once_with(cube, "dummy_fnme")

    def test_no_saver_found(self):
        """Unrecognised saver specified"""
        msg = "Cannot save; no saver can be found associated with " '"dummy_saver"'
        with self.assertRaisesRegex(ValueError, msg):
            with mock.patch("warnings.warn") as mock_warn:
                ants.save(self.cube, "dummy_fnme", saver="dummy_saver")
            self.assertEqual(mock_warn.call_count, 1)
            calls = [ants_save_call()]
            mock_warn.assert_has_calls(calls)

    def test_fill_value_set_correctly_in_netcdf_file(self):
        fh = tempfile.NamedTemporaryFile(suffix=".nc")
        today = datetime.datetime(1, 1, 1)
        with mock.patch("ants.fileformats.datetime") as dt_patch:
            dt_patch.today.return_value = today
            with mock.patch("sys.argv", new=["program", "arg1", "arg2"]):
                cube = stock.lat_lon_cube()
                cube.attributes["STASH"] = "m01s01i001"
                with mock.patch("warnings.warn") as mock_warn:
                    ants.save(cube, fh.name, fill_value=-99)
                self.assertEqual(mock_warn.call_count, 2)
                calls = [ants_save_call(), netcdf_save_call()]
                mock_warn.assert_has_calls(calls)
        self.assertCDL(fh.name, ("fileformats", "fill_value_set_netcdf.cdl"))

    def test_fill_value_handled_correctly_by_ancil_save(self):
        # test to ensure passing a fill_value to ants save and
        # specifying "ancil" saver does not throw an error
        anc_fh = tempfile.NamedTemporaryFile()
        self.cdl_filename_to_delete = f"{anc_fh.name}.nc"
        today = datetime.datetime(1, 1, 1)
        with mock.patch("ants.fileformats.datetime") as dt_patch:
            dt_patch.today.return_value = today
            with mock.patch("sys.argv", new=["program", "arg1", "arg2"]):
                cube = stock.lat_lon_cube()
                cube.attributes["STASH"] = "m01s01i001"
                with mock.patch("warnings.warn") as mock_warn:
                    ants.save(cube, anc_fh.name, fill_value=-99, saver="ancil")
                self.assertEqual(mock_warn.call_count, 3)
                calls = [ants_save_call(), netcdf_save_call(), ancil_save_call()]
                mock_warn.assert_has_calls(calls)
        self.assertCDL(f"{anc_fh.name}.nc", ("fileformats", "ancillary_file.cdl"))


class TestLoadCube(ants.tests.TestCase):
    def test_circular_derive(self):
        # Ensure that we can derive global field status on load (circular)
        cube = ants.tests.stock.geodetic((2, 2))
        cube.coord(axis="x").circular = False
        fh = tempfile.NamedTemporaryFile(suffix=".nc")
        with mock.patch("warnings.warn") as mock_warn:
            ants.save(cube, fh.name)
        self.assertEqual(mock_warn.call_count, 2)
        calls = [ants_save_call(), netcdf_save_call()]
        mock_warn.assert_has_calls(calls)
        result = ants.load_cube(fh.name)
        self.assertTrue(result.coord(axis="x").circular)

    def test_no_horizontal_coords(self):
        # Ensure that in the case where we have no horizontal coords, that
        # we are still able to load a cube.
        cube = iris.cube.Cube([1])
        with mock.patch("iris.io.load_files", return_value=[cube]):
            ants.load_cube("stub_filename")

    @ants.tests.skip_gdal
    def test_raster_load(self):
        cube = ants.tests.stock.geodetic((6, 3))
        res_cube = ants.load_cube(ants.tests.get_data_path("global_geodetic.bil"))
        self.assertEqual(res_cube, cube)

    def test_namelist_load(self):
        fh = tempfile.NamedTemporaryFile()
        data = b"""&GRID
 POINTS_LAMBDA_TARG=2,POINTS_PHI_TARG=2
/
"""
        fh.write(data)
        fh.seek(0)
        CAPGridRegular_dummy = mock.MagicMock()
        CAPGridRegular_dummy.get_cube.return_value = ants.tests.stock.geodetic((2, 2))
        patch = mock.patch(
            "ants.fileformats.namelist.CAPGridRegular",
            return_value=CAPGridRegular_dummy,
        )
        with patch as npatch:
            ants.load_cube(fh.name)
        self.assertTrue(npatch.called)

    def test_forecast_period_removal(self):
        cube = ants.tests.stock.geodetic((2, 2))
        cube.add_aux_coord(iris.coords.AuxCoord(0, long_name="forecast_period"))
        fh = tempfile.NamedTemporaryFile(suffix=".nc")
        with mock.patch("warnings.warn") as mock_warn:
            ants.save(cube, fh.name)
        self.assertEqual(mock_warn.call_count, 2)
        calls = [ants_save_call(), netcdf_save_call()]
        mock_warn.assert_has_calls(calls)
        result = ants.load_cube(fh.name)
        result_coords = [c.name() for c in result.coords()]
        self.assertNotIn("forecast_period", result_coords)

    def test_forecast_reference_time_removal(self):
        cube = ants.tests.stock.geodetic((2, 2))
        cube.add_aux_coord(iris.coords.AuxCoord(0, long_name="forecast_reference_time"))
        fh = tempfile.NamedTemporaryFile(suffix=".nc")
        with mock.patch("warnings.warn") as mock_warn:
            ants.save(cube, fh.name)
        self.assertEqual(mock_warn.call_count, 2)
        calls = [ants_save_call(), netcdf_save_call()]
        mock_warn.assert_has_calls(calls)
        result = ants.load_cube(fh.name)
        result_coords = [c.name() for c in result.coords()]
        self.assertNotIn("forecast_reference_time", result_coords)


class TestPriorities(ants.tests.TestCase):
    def test_grib(self):
        grib_format_specs = [
            spec
            for spec in iris.fileformats.FORMAT_AGENT._format_specs
            if spec.name == "GRIB"
        ]
        if len(grib_format_specs) != 1:
            self.fail(f"{len(grib_format_specs)} GRIB format specs found, expected 1.")
        grib_format_spec = grib_format_specs[0]

        gdal_format_specs = [
            spec
            for spec in iris.fileformats.FORMAT_AGENT._format_specs
            if spec.name == "gdal"
        ]
        if len(gdal_format_specs) != 1:
            self.fail(f"{len(gdal_format_specs)} GDAL format specs found, expected 1.")
        gdal_format_spec = gdal_format_specs[0]

        self.assertGreater(grib_format_spec.priority, gdal_format_spec.priority)


if __name__ == "__main__":
    ants.tests.main()
