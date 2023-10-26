# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import tempfile
import unittest.mock as mock

import ants
import ants.io.save as save
import ants.tests
import ants.tests.stock
import iris


class TestLoadCube(ants.tests.TestCase):
    def test_circular_derive(self):
        # Ensure that we can derive global field status on load (circular)
        cube = ants.tests.stock.geodetic((2, 2))
        cube.coord(axis="x").circular = False
        fh = tempfile.NamedTemporaryFile(suffix=".nc")
        save.netcdf(cube, fh.name)
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
        save.netcdf(cube, fh.name)
        result = ants.load_cube(fh.name)
        result_coords = [c.name() for c in result.coords()]
        self.assertNotIn("forecast_period", result_coords)

    def test_forecast_reference_time_removal(self):
        cube = ants.tests.stock.geodetic((2, 2))
        cube.add_aux_coord(iris.coords.AuxCoord(0, long_name="forecast_reference_time"))
        fh = tempfile.NamedTemporaryFile(suffix=".nc")
        save.netcdf(cube, fh.name)
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
