# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import tempfile
import unittest
from unittest import mock

import ants.io.save as save
import ants.tests
import ants.tests.stock


class TestAll(ants.tests.TestCase):
    def test_netcdf_classic_int64_coordinates(self):
        # NETCDF_CLASSIC does not support 64bit integer coordinates, ensure
        # that we can handle them.
        cube = ants.tests.stock.geodetic((2, 2))
        coord = cube.coord(axis="x")
        coord.points = coord.points.astype("int64")
        fh = tempfile.NamedTemporaryFile()
        self.assertIsNone(save.netcdf(cube, fh.name))

    def test_netcdf_classic_int64_data(self):
        # NETCDF_CLASSIC does not support 64bit integer data, ensure
        # that we can handle them (via iris).
        cube = ants.tests.stock.geodetic((2, 2))
        cube.data = cube.data.astype("int64")
        fh = tempfile.NamedTemporaryFile()
        self.assertIsNone(save.netcdf(cube, fh.name))

    def test_netcdf_classic_int64_data_lazy_with_ants_override(self):
        # Test that the ANTS workaround for NETCDF4_CLASSIC saving with
        # lazy int64 data enables the save to happen.
        # See https://code.metoffice.gov.uk/trac/ancil/ticket/1527
        cube = ants.tests.stock.geodetic((2, 2))
        cube = ants.utils.cube.defer_cube(cube)
        cube.data = cube.core_data().astype("int64")
        fh = tempfile.NamedTemporaryFile()
        self.assertIsNone(save.netcdf(cube, fh.name))

    @unittest.expectedFailure
    @mock.patch("ants.io.save._iris_netcdf4_classic_workaround")
    def test_netcdf_classic_int64_data_lazy_without_ants_override(self, *args):
        # Test that the ANTS workaround for NETCDF4_CLASSIC saving with
        # lazy int64 data is needed. If this test passes, the workaround
        # can be removed. This should happen when ANTS is upgraded to Iris
        # 3.2.
        # See https://code.metoffice.gov.uk/trac/ancil/ticket/1527
        cube = ants.tests.stock.geodetic((2, 2))
        cube = ants.utils.cube.defer_cube(cube)
        cube.data = cube.core_data().astype("int64")
        fh = tempfile.NamedTemporaryFile()
        self.assertIsNone(save.netcdf(cube, fh.name))


if __name__ == "__main__":
    ants.tests.main()
