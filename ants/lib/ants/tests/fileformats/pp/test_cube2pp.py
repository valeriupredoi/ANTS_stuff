# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import unittest
import unittest.mock as mock

import ants
import ants.fileformats.pp as pp
import ants.tests
import cf_units
import iris
import numpy as np


class TestYearZero(ants.tests.TestCase):
    def test_supplemented_error_msg(self):
        # We supplement date handling errors by pointing the user to the
        # appendix for further guidance.
        cube = ants.tests.stock.geodetic((2, 2, 2))
        points = [92, 274]
        bounds = [[0, 183], [183, 366]]
        tunit = cf_units.Unit("days since 0000-01-01 00:00:00", calendar="gregorian")
        time_coord = iris.coords.DimCoord(
            points, bounds=bounds, standard_name="time", units=tunit
        )
        cube.add_dim_coord(time_coord, 0)
        msg = "zero not allowed.*\nEnsure that dates are representative.*"
        with self.assertRaisesRegex(ValueError, msg):
            pp.cube2pp(cube)


class TestZonalMean_ants_no_override(ants.tests.TestCase):
    """
    These tests assume correct iris behaviour, which may or may not be valid.
    Check TestZonalMean_iris_failure below for details of the current state of
    iris responses.

    """

    def setUp(self):
        self.cube = ants.tests.stock.geodetic((2, 1))
        self.field = mock.Mock(
            "field",
            lbnpt=1,
            lbuser=[2],
            bdy=180.0,
            bdx=360.0,
            lbproc=64.0,
            bplat=90,
            bplon=0,
            lblev=9999,
            lbcode=mock.Mock("lbcode", ix=0),
        )

    def run_cube2pp(self, field):
        with mock.patch("iris.fileformats.pp.as_fields", return_value=[field]):
            return pp.cube2pp(self.cube)

    def test_lbnpt(self):
        fields = self.run_cube2pp(self.field)
        self.assertEqual(fields[0].lbnpt, 1)

    def test_lbproc(self):
        fields = self.run_cube2pp(self.field)
        self.assertEqual(fields[0].lbproc, 64)

    def test_bdx(self):
        fields = self.run_cube2pp(self.field)
        self.assertEqual(fields[0].bdx, 360.0)


class TestZonalMean_ants_override(TestZonalMean_ants_no_override):
    """
    These tests ensure that our workarounds behave as expected for current
    iris treatment of zonal mean data.

    """

    def setUp(self):
        self.cube = ants.tests.stock.geodetic((2, 1))
        self.field = mock.Mock(
            "field",
            lbnpt=0,
            bdy=1,
            bdx=0.0,
            lbuser=[2],
            lbproc=64,
            bplat=90,
            bplon=0,
            lblev=9999,
            lbcode=mock.Mock("lbcode", ix=0),
        )

    def test_bdx_with_bdy_missing(self):
        field = mock.Mock(
            "field",
            lbnpt=0,
            bdy=pp.RMDI,
            bdx=0.0,
            lbuser=[2],
            lbproc=64,
            bplat=90,
            bplon=0,
            lblev=9999,
            lbcode=mock.Mock("lbcode", ix=0),
        )
        fields = self.run_cube2pp(field)
        self.assertEqual(fields[0].bdx, 0.0)


class TestZonalMean_iris_behaviour(TestZonalMean_ants_override):
    """
    These tests ensure that our workarounds behave as expected for current
    iris treatment of zonal mean data.

    """

    def run_cube2pp(self, field):
        with mock.patch("ants.fileformats.pp._iris_workarounds"):
            return super().run_cube2pp(field)

    # Both test workarounds in cube2pp.py can be removed when iris fixed is
    # available.  Iris issue: https://github.com/SciTools/iris/issues/2054
    @unittest.expectedFailure
    def test_lbnpt(self):
        return super().test_lbnpt()

    @unittest.expectedFailure
    def test_bdx(self):
        return super().test_bdx()


class Test_unrotated_lam_behaviour(ants.tests.TestCase):
    """Check cube2pp behaviour for unrotated pole LAMs."""

    def setUp(self):
        self.cube = ants.tests.stock.geodetic((2, 2), xlim=(0, 90))

    def check_values(self):
        fields = pp.cube2pp(self.cube)
        self.assertEqual(180, fields[0].bplon)

    def test_ants_override(self):
        self.check_values()

    @unittest.expectedFailure
    def test_ants_no_override(self):
        # If this succeeds then we can remove this test class and
        # corresponding workaround in ants.fileformats.pp._iris_workarounds
        with mock.patch("ants.fileformats.pp._iris_workarounds"):
            self.check_values()


class TestSetLbuser0(ants.tests.TestCase):
    """
    Normally lbuser0 isn't set until actually writting to disk by iris.  We
    manually set this so that we can correctly declare boolean fields but also
    so that F03 ancillaries derived from iris pp fields can simply inherrit the
    data type metadata.

    """

    def test_int_in_int_out(self):
        cube = ants.tests.stock.geodetic(data=np.zeros((2, 2), dtype="int"))
        field = list(pp.cube2pp(cube))[0]
        self.assertEqual(field.lbuser[0], 2)

    def test_uint_in_int_out(self):
        cube = ants.tests.stock.geodetic(data=np.zeros((2, 2), dtype="uint8"))
        field = list(pp.cube2pp(cube))[0]
        self.assertEqual(field.lbuser[0], 2)

    def test_int_in_logical_out(self):
        # Where valid_range or valid_min & valid_max are provided and
        # correspond to [0, 1], this denotes a logical field.
        cube = ants.tests.stock.geodetic(data=np.zeros((2, 2), dtype="int"))
        cube.attributes = {"valid_range": [0, 1]}
        field = list(pp.cube2pp(cube))[0]
        self.assertEqual(field.lbuser[0], 3)

    def test_uint_in_logical_out(self):
        # Where valid_range or valid_min & valid_max are provided and
        # correspond to [0, 1], this denotes a logical field.
        cube = ants.tests.stock.geodetic(data=np.zeros((2, 2), dtype="uint8"))
        cube.attributes = {"valid_range": [0, 1]}
        field = list(pp.cube2pp(cube))[0]
        self.assertEqual(field.lbuser[0], 3)

    def test_bool_in_logical_out(self):
        cube = ants.tests.stock.geodetic(data=np.zeros((2, 2), dtype="bool"))
        field = list(pp.cube2pp(cube))[0]
        self.assertEqual(field.lbuser[0], 3)

    def test_float_in_float_out(self):
        cube = ants.tests.stock.geodetic(data=np.zeros((2, 2), dtype="float"))
        field = list(pp.cube2pp(cube))[0]
        self.assertEqual(field.lbuser[0], 1)


class TestValidity(ants.tests.TestCase):
    def test_overspecified(self):
        cube = ants.tests.stock.geodetic((2, 2))
        cube.attributes = {
            "valid_range": [0, 1],
            "valid_max": 2,
            "valid_min": 1,
        }
        msg = "Cube attribute valid range is overspecified"
        with self.assertRaisesRegex(ValueError, msg):
            pp.cube2pp(cube)


class TestDefinitionSouthNorth(ants.tests.TestCase):
    def setUp(self):
        self.cube = ants.tests.stock.geodetic((2, 2), ylim=(90, -90))

    def test_original_unchanged(self):
        expected = self.cube.copy()
        pp.cube2pp(self.cube)
        self.assertEqual(expected, self.cube)

    def test_south_north(self):
        # Needs return value so remainder of pp.cube2pp doesn't throw
        # exception.
        with mock.patch(
            "ants.utils.cube.reverse_coordinate", return_value=self.cube
        ) as patch:
            pp.cube2pp(self.cube)
        patch.assert_called_once_with(self.cube, self.cube.coord(axis="Y"))


class TestTimeCoordinateHandling(ants.tests.TestCase):
    def test_multiple_time_coordinates(self):
        # We supplement date handling errors by pointing the user to the
        # appendix for further guidance.
        cube = ants.tests.stock.geodetic((2, 2, 2))
        points = [92, 274]
        bounds = [[0, 183], [183, 366]]
        tunit = cf_units.Unit("days since 0000-01-01 00:00:00", calendar="gregorian")
        time_coord = iris.coords.DimCoord(
            points, bounds=bounds, standard_name="time", units=tunit
        )
        extra_time_coord = iris.coords.DimCoord(
            points,
            bounds=bounds,
            standard_name="time",
            long_name="extra_time_coord",
            units=tunit,
        )
        cube.add_aux_coord(time_coord, 0)
        cube.add_aux_coord(extra_time_coord, 0)
        msg = "More than one time based coordinate"
        with self.assertRaisesRegex(RuntimeError, msg):
            pp.cube2pp(cube)

    def test_gregorian_calendar(self):
        # Check behaviour for gregorian calendar.
        cube = ants.tests.stock.geodetic((2, 2, 2))
        points = [92, 274]
        bounds = [[0, 183], [183, 366]]
        tunit = cf_units.Unit("days since 0001-01-01 00:00:00", calendar="gregorian")
        time_coord = iris.coords.DimCoord(
            points, bounds=bounds, standard_name="time", units=tunit
        )
        cube.add_aux_coord(time_coord, 0)
        result = pp.cube2pp(cube)
        self.assertEqual(1, result[0].lbtim)

    def test_proleptic_gregorian_calendar(self):
        # Check that proleptic gregorian calendar is being correctly converted.
        cube = ants.tests.stock.geodetic((2, 2, 2))
        points = [92, 274]
        bounds = [[0, 183], [183, 366]]
        tunit = cf_units.Unit(
            "days since 0001-01-01 00:00:00", calendar="proleptic_gregorian"
        )
        time_coord = iris.coords.DimCoord(
            points, bounds=bounds, standard_name="time", units=tunit
        )
        cube.add_aux_coord(time_coord, 0)
        result = pp.cube2pp(cube)
        self.assertEqual(1, result[0].lbtim)


class TestSetLblevSurfaceField(ants.tests.TestCase):
    """
    Check that the lblev value for surface fields is set to 9999
    rather than 0.

    """

    def setUp(self):
        self.cube = ants.tests.stock.geodetic((2, 2))

    def check_values(self):
        fields = pp.cube2pp(self.cube)
        self.assertEqual(9999, fields[0].lblev)

    def test_ants_override(self):
        self.check_values()

    @unittest.expectedFailure
    def test_ants_no_override(self):
        # If this succeeds then we can remove this test class and
        # corresponding workaround in
        # ants.fileformats.pp._set_lblev_surface_field. Iris issue:
        # https://github.com/SciTools/iris/issues/3820
        with mock.patch("ants.fileformats.pp._set_lblev_surface_field"):
            self.check_values()


class TestSetLblevDepthField(ants.tests.TestCase):
    """
    Check that the lblev value for depth fields is set to model level number
    rather than 0 (which is set to 9999).

    """

    def setUp(self):
        self.cube = ants.tests.stock.geodetic((2, 2, 2))
        points = [0.5, 2.0]
        bounds = [[0.0, 1.0], [183, 366]]
        dunit = cf_units.Unit("m")
        depth_coord = iris.coords.DimCoord(
            points,
            bounds=bounds,
            standard_name="depth",
            units=dunit,
            attributes={"positive": "down"},
        )
        self.cube.add_dim_coord(depth_coord, 0)

    def assert_lblevs(self):
        fields = pp.cube2pp(self.cube)
        self.assertEqual(1, fields[0].lblev)
        self.assertEqual(2, fields[1].lblev)

    def test_ants_override(self):
        self.assert_lblevs()

    @unittest.expectedFailure
    def test_ants_no_override(self):
        # If this succeeds then we can remove this test class and
        # corresponding workaround in
        # ants.fileformats.pp._set_lblev_depths. Iris issue:
        # https://github.com/SciTools/iris/issues/4082.
        with mock.patch("ants.fileformats.pp._set_lblev_depths"):
            self.assert_lblevs()


if __name__ == "__main__":
    ants.tests.main()
