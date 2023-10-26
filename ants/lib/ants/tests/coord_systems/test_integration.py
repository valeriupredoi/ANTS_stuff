# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.coord_systems as coord_systems
import ants.tests
import iris
import numpy as np
from iris.coord_systems import OSGB, GeogCS


class Test_ants_coordinate_systems(ants.tests.TestCase):
    def test_wgs84_as_umsphere(self):
        # Coordinate systems considered equivalent.
        wgs84 = coord_systems.WGS84_GEODETIC.crs
        wgs84_v2 = GeogCS(semi_major_axis=6378137.0, inverse_flattening=298.257223563)
        res_crs = wgs84.as_ants_crs()
        self.assertEqual(res_crs, coord_systems.UM_SPHERE.crs)

        res_crs = wgs84_v2.as_ants_crs()
        self.assertEqual(res_crs, coord_systems.UM_SPHERE.crs)

    def test_osgb_as_transverse_mercator(self):
        # Ensure that OSGB returns the more general Transverse Mercator in
        # order to remove projection limits in its usage.
        osgb = OSGB()
        res_crs = osgb.as_ants_crs()
        self.assertEqual(res_crs, coord_systems.OSGB.crs)
        self.assertEqual(type(res_crs), iris.coord_systems.TransverseMercator)

    def test_identical_return(self):
        # Coordinate systems with no equivalent
        crs = GeogCS(6000000)
        res_crs = crs.as_ants_crs()
        self.assertIs(res_crs, crs)

    def test_rotated_pole_treatment_as_unrotated(self):
        crs = iris.coord_systems.RotatedGeogCS(
            90, 180, ellipsoid=coord_systems.UM_SPHERE.crs
        )
        res_crs = crs.as_ants_crs()
        self.assertIsInstance(res_crs, coord_systems.UM_SPHERE.crs.__class__)

    def test_rotated_pole_treatment_unchanged(self):
        crs = iris.coord_systems.RotatedGeogCS(
            30, 180, ellipsoid=coord_systems.UM_SPHERE.crs
        )
        res_crs = crs.as_ants_crs()
        self.assertIs(res_crs, crs)


class Test_crs_float_equality(ants.tests.TestCase):
    # If/when iris handles coordinate system floating comparison correctly,
    # we can remove this workaround of overriding the crs equality special
    # method.
    # These tests are the basis for the chosen tolerance when comparison
    # coordinate systems.  Note that 32bit floats end up on coordinate
    # system attributes via pp loading.

    def test_derived_quantities_eq(self):
        # Ensure that derived quantities are equal.
        crs1 = GeogCS(semi_major_axis=6377563.396, semi_minor_axis=6356256.909)
        crs2 = GeogCS(semi_major_axis=6377563.396, inverse_flattening=299.324)
        self.assertEqual(crs1, crs2)

    def test_32bit_derived_quantities_eq(self):
        # Ensure that 32-bit coerced quantities don't fail equality.
        crs1 = GeogCS(semi_major_axis=6377563.396, semi_minor_axis=6356256.909)
        crs2 = GeogCS(
            semi_major_axis=np.float32(6377563.396),
            inverse_flattening=np.float32(299.324),
        )
        self.assertEqual(crs1, crs2)

    def test_derived_quantities_neq(self):
        # Ensure that derived quantities are not equal.
        crs1 = GeogCS(semi_major_axis=6377563.396, semi_minor_axis=6356256.909)
        crs2 = GeogCS(semi_major_axis=6377563.396, inverse_flattening=299.33)
        self.assertFalse(crs1 == crs2)

    def test_straight_eq(self):
        # No derived quantities.
        crs1 = GeogCS(6377560)
        crs2 = GeogCS(6377570)
        self.assertEqual(crs1, crs2)

    def test_straight_neq(self):
        # No derived quantities.
        crs1 = GeogCS(6377500)
        crs2 = GeogCS(6377600)
        self.assertFalse(crs1 == crs2)

    def test_straight_32bit_eq(self):
        # No derived quantities.
        crs1 = GeogCS(6377560)
        crs2 = GeogCS(np.float32(6377570))
        self.assertEqual(crs1, crs2)

    def test_straight_32bit_neq(self):
        # No derived quantities.
        crs1 = GeogCS(6377500)
        crs2 = GeogCS(np.float32(6377600))
        self.assertFalse(crs1 == crs2)


if __name__ == "__main__":
    ants.tests.main()
