# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
import numpy as np
from ants.regrid._ugrid import LATITUDE_DIM, LONGITUDE_DIM
from ants.regrid._ugrid import _build_grid as build_grid


@ants.tests.skip_esmf
class TestAll(ants.tests.TestCase):
    def test_type(self):
        import ESMF

        expected = ESMF.Grid

        actual = build_grid(ants.tests.stock.geodetic((2, 2)))

        self.assertIsInstance(actual, expected)

    def test_grid_from_multidimensional_source(self):
        # Check get the same shape of output for 2D and 3D cube as long as
        # they have the same horizontal coordinates:
        expected = build_grid(ants.tests.stock.geodetic((5, 6)))

        actual = build_grid(ants.tests.stock.simple_3d_time_varying())

        self.assertArrayEqual(expected.max_index, actual.max_index)

    def test_grid_from_masked_multidimensional_source(self):
        # Check get the same shape of output for 2D and 3D cube as long as
        # they have the same horizontal coordinates:
        expected = build_grid(ants.tests.stock.geodetic((5, 6)))

        masked_cube = ants.tests.stock.simple_3d_time_varying()
        masked_cube.data = np.ma.masked_array(masked_cube.data)
        masked_cube.data[..., 0] = np.ma.masked
        actual = build_grid(masked_cube)

        self.assertArrayEqual(expected.max_index, actual.max_index)

    def test_masked_with_rotated_pole_coordinates(self):
        shape = (3, 3, 3)
        cube = ants.tests.stock.geodetic(shape, north_pole_lat=80, north_pole_lon=80)
        mask = np.zeros(shape, dtype=np.bool)
        mask[0, 0, 0] = True
        cube.data = np.ma.masked_array(cube.data, mask=mask)

        self.assertIsNotNone(build_grid(cube))


@ants.tests.skip_esmf
class TestGridCoordinates(ants.tests.TestCase):
    # Tests for multiple longitudes (i.e. not zonal mean data).  Note that
    # returned coordinates are the coordinate value as a complete array not as
    # an axis.  So for a 4x3 grid, there's a 4x3 array of latitude centre
    # values and a 4x3 array of longitude centre values.  For the corners,
    # there's one extra value in each dimension (i.e. 5x4 corner locations for
    # 4x3 grid).

    # 4x3 example:
    #
    # 1 --- 2 --- 3 --- 4 --- 5
    # |     |     |     |     |
    # |  A  |  B  |  C  |  D  |
    # |     |     |     |     |
    # 6 --- 7 --- 8 --- 9 --- 10
    # |     |     |     |     |
    # |  E  |  F  |  G  |  H  |
    # |     |     |     |     |
    # 11 ---12 ---13 ---14 ---15
    # |     |     |     |     |
    # |  I  |  J  |  K  |  L  |
    # |     |     |     |     |
    # 16 ---17 ---18 ---19 ---20
    #
    # 4x3 cell centre locations, labelled A through L.
    #
    # 5x4 cell corner locations, labelled 1 through 20.
    def setUp(self):
        self.columns = 4
        self.rows = 3

    def test_grid_corner_longitudes(self):
        import ESMF

        longitudes = np.linspace(-180.0, 180.0, self.columns + 1)
        expected = np.tile(longitudes, (self.rows + 1, 1))

        grid = build_grid(ants.tests.stock.geodetic((self.rows, self.columns)))
        actual = grid.get_coords(
            LONGITUDE_DIM, staggerloc=ESMF.api.constants.StaggerLoc.CORNER
        )

        self.assertArrayAlmostEqual(expected, actual)

    def test_grid_centre_longitudes(self):
        import ESMF

        longitudes = np.linspace(-135, 135, self.columns)
        expected = np.tile(longitudes, (self.rows, 1))

        grid = build_grid(ants.tests.stock.geodetic((self.rows, self.columns)))
        actual = grid.get_coords(
            LONGITUDE_DIM, staggerloc=ESMF.api.constants.StaggerLoc.CENTER
        )
        self.assertArrayAlmostEqual(expected, actual)

    def test_grid_corner_latitudes(self):
        import ESMF

        latitudes = np.linspace(-90.0, 90.0, self.rows + 1)
        expected = np.tile(latitudes, (self.columns + 1, 1)).T

        grid = build_grid(ants.tests.stock.geodetic((self.rows, self.columns)))
        actual = grid.get_coords(
            LATITUDE_DIM, staggerloc=ESMF.api.constants.StaggerLoc.CORNER
        )

        self.assertArrayAlmostEqual(expected, actual)

    def test_grid_centre_latitudes(self):
        import ESMF

        latitudes = np.linspace(-60, 60, self.rows)
        expected = np.tile(latitudes, (self.columns, 1)).T

        grid = build_grid(ants.tests.stock.geodetic((self.rows, self.columns)))
        actual = grid.get_coords(
            LATITUDE_DIM, staggerloc=ESMF.api.constants.StaggerLoc.CENTER
        )
        self.assertArrayAlmostEqual(expected, actual)


@ants.tests.skip_esmf
class TestZonalMeanCoordinates(ants.tests.TestCase):
    # Note that returned coordinates are the coordinate value as a complete
    # array not as an axis.  So for a 1x3 grid, there's a 1x3 array of
    # latitude centre values and a 1x3 array of longitude centre values.  For
    # the corners, there's one extra value in each dimension (i.e. 2x4 corner
    # locations for 1x3 grid).
    #
    # 1 --------------------- 2
    # |                       |
    # |           A           |
    # |                       |
    # 3 --------------------- 4
    # |                       |
    # |           B           |
    # |                       |
    # 5 --------------------- 6
    # |                       |
    # |           C           |
    # |                       |
    # 7 --------------------- 8
    #
    # 1x3 cell centre locations, labelled A through C.
    #
    # 2x4 cell corner locations, labelled 1 through 8.
    def setUp(self):
        self.columns = 1
        self.rows = 3
        # The fix converts single zonal mean longitude into 12 longitudes:
        self.zonal_mean_fix_columns = 12

    def test_grid_corner_longitudes_from_zonal_mean(self):
        import ESMF

        longitudes = np.linspace(-180.0, 180.0, self.zonal_mean_fix_columns + 1)
        expected = np.tile(longitudes, (self.rows + 1, 1))

        grid = build_grid(ants.tests.stock.geodetic((self.rows, self.columns)))
        actual = grid.get_coords(
            LONGITUDE_DIM, staggerloc=ESMF.api.constants.StaggerLoc.CORNER
        )

        self.assertArrayAlmostEqual(expected, actual)

    def test_grid_centre_longitudes_from_zonal_mean(self):
        import ESMF

        longitudes = np.linspace(-165.0, 165.0, self.zonal_mean_fix_columns)
        expected = np.tile(longitudes, (self.rows, 1))

        grid = build_grid(ants.tests.stock.geodetic((self.rows, self.columns)))
        actual = grid.get_coords(
            LONGITUDE_DIM, staggerloc=ESMF.api.constants.StaggerLoc.CENTER
        )
        self.assertArrayAlmostEqual(expected, actual)

    def test_grid_corner_latitudes_from_zonal_mean(self):
        import ESMF

        latitudes = np.linspace(-90.0, 90.0, self.rows + 1)
        expected = np.tile(latitudes, (self.zonal_mean_fix_columns + 1, 1)).T

        grid = build_grid(ants.tests.stock.geodetic((self.rows, self.columns)))
        actual = grid.get_coords(
            LATITUDE_DIM, staggerloc=ESMF.api.constants.StaggerLoc.CORNER
        )

        self.assertArrayAlmostEqual(expected, actual)

    def test_grid_centre_latitudes_from_zonal_mean(self):
        import ESMF

        latitudes = np.linspace(-60, 60, self.rows)
        expected = np.tile(latitudes, (self.zonal_mean_fix_columns, 1)).T

        grid = build_grid(ants.tests.stock.geodetic((self.rows, self.columns)))
        actual = grid.get_coords(
            LATITUDE_DIM, staggerloc=ESMF.api.constants.StaggerLoc.CENTER
        )
        self.assertArrayAlmostEqual(expected, actual)

    def test_center_staggerloc_mask(self):
        expected = np.array([[1, 0], [0, 0]], dtype=np.int32)

        source = ants.tests.stock.geodetic((2, 2))
        source.data = np.ma.masked_array(source.data)
        source.data[0, 0] = np.ma.masked
        grid = build_grid(source)
        # Note: always get all 4 staggerlocs (centre, edge1, edge2 and
        # corner) back from ESMPy Grid:
        # http://www.earthsystemmodeling.org/esmf_releases/public/ESMF_7_1_0r/esmpy_doc/html/StaggerLoc.html#ESMF.api.constants.StaggerLoc
        # but we only want to mask centre.
        actual = grid.mask[0]

        self.assertArrayEqual(expected, actual)

    def test_other_staggerlocs_not_masked(self):
        expected = [None, None, None]

        source = ants.tests.stock.geodetic((2, 2))
        source.data = np.ma.masked_array(source.data)
        source.data[0, 0] = np.ma.masked
        grid = build_grid(source)
        # Note: always get all 4 staggerlocs (centre, edge1, edge2 and
        # corner) back from ESMPy Grid:
        # http://www.earthsystemmodeling.org/esmf_releases/public/ESMF_7_1_0r/esmpy_doc/html/StaggerLoc.html#ESMF.api.constants.StaggerLoc
        # but we only want to mask centre.
        actual = grid.mask[1:]

        self.assertEqual(expected, actual)


if __name__ == "__main__":
    ants.tests.main()
