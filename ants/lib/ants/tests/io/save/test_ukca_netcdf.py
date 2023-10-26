# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import unittest.mock as mock

import ants.io.save as save
import ants.tests
import iris.tests.stock as stock
import numpy as np


class TestAll(ants.tests.TestCase):
    def setUp(self):
        self.filename = "my_filename"
        self.cubes = [ants.tests.stock.geodetic((3, 2))]
        self.cubes[0].data = np.ma.array(self.cubes[0].data)
        self.cubes[0].data[:, 0] = np.ma.masked
        self.cubes[0].attributes["emission_type"] = "2"

    def test_default_args(self):
        with mock.patch("iris.save") as mock_save:
            save.ukca_netcdf(self.cubes, self.filename)
        mock_save.assert_called_once_with(
            self.cubes,
            f"{self.filename}.nc",
            saver="nc",
            netcdf_format="NETCDF4_CLASSIC",
            local_keys=[
                "tracer_name",
                "vertical_scaling",
                "lowest_level",
                "highest_level",
                "hourly_scaling",
            ],
            unlimited_dimensions=None,
            zlib=True,
            complevel=4,
            fill_value=None,
        )

    def test_ukca_conventions_applied_on_realised_data(self):
        with mock.patch("iris.save"):
            save.ukca_netcdf(self.cubes, self.filename)
        self.assertEqual(self.cubes[0].data.dtype, np.int32)
        self.assertEqual(self.cubes[0].attributes["update_type"].dtype, np.int32)
        self.assertNotIn("emission_type", self.cubes[0].attributes)
        target_data = np.arange(6).reshape((3, 2))
        # Mask replaced with unmasked numpy array with 0.
        self.assertFalse(np.ma.isMaskedArray(self.cubes[0].data))
        target_data[:, 0] = 0
        self.assertArrayEqual(self.cubes[0].data, target_data)

    def test_ukca_conventions_applied_on_lazy_data(self):
        # Need to defer the cube so that the data is lazy.
        cube = ants.utils.cube.defer_cube(self.cubes[0])
        with mock.patch("iris.save"):
            save.ukca_netcdf(cube, self.filename)
        # Check that the cube.data is still lazy.
        self.assertTrue(cube.has_lazy_data())
        self.assertEqual(cube.data.dtype, np.int32)
        self.assertEqual(cube.attributes["update_type"].dtype, np.int32)
        self.assertNotIn("emission_type", cube.attributes)
        target_data = np.arange(6).reshape((3, 2))
        # Mask replaced with unmasked numpy array with 0.
        self.assertFalse(np.ma.isMaskedArray(cube.data))
        target_data[:, 0] = 0
        self.assertArrayEqual(cube.data, target_data)


class TestSave(ants.tests.TestCase):
    def setUp(self):
        self.cube = stock.lat_lon_cube()

    def test_ukca_netcdf_save(self):
        ukcapatch = mock.patch("ants.io.save._ukca_conventions")
        cfpatch = mock.patch("ants.io.save.netcdf")
        with cfpatch, ukcapatch as nc_patch:
            save.ukca_netcdf(self.cube, "dummy_fnme")
        nc_patch.assert_called_once()


if __name__ == "__main__":
    ants.tests.main()
