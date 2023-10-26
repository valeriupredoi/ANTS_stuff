# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import unittest.mock as mock

import ants.tests
from ants.regrid.rectilinear import TwoStage


class Testall(ants.tests.TestCase):
    def setUp(self):
        patch = mock.patch("ants.regrid.rectilinear._AreaWeightedRegridder")
        self.area_patch = patch.start()
        self.addCleanup(patch.stop)

        patch = mock.patch("ants.regrid.rectilinear.Linear.regridder")
        self.linear_patch = patch.start()
        self.addCleanup(patch.stop)

        patch = mock.patch("ants.regrid.rectilinear._gen_regular_target")
        self.regtar_patch = patch.start()
        self.addCleanup(patch.stop)

    def test_same_crs_regrid(self):
        src_cube = ants.tests.stock.geodetic((4, 4))
        tgt_cube = ants.tests.stock.geodetic((4, 4))
        regridder = TwoStage()
        re = regridder.regridder(src_cube, tgt_cube)
        re(src_cube)

        self.assertTrue(self.area_patch.called)
        self.assertFalse(self.regtar_patch.called)
        self.assertFalse(self.linear_patch.called)

    def test_diff_crs_regrid(self):
        src_cube = ants.tests.stock.geodetic((4, 4))
        tgt_cube = ants.tests.stock.osgb((4, 4))
        regridder = TwoStage()
        re = regridder.regridder(src_cube, tgt_cube)
        re(src_cube)

        self.assertTrue(self.area_patch.called)
        self.assertTrue(self.regtar_patch.called)
        self.assertTrue(self.linear_patch.called)

    def test_homogenised_crs(self):
        src_cube = ants.tests.stock.geodetic((4, 4))
        crs = ants.coord_systems.WGS84_GEODETIC.crs
        src_cube.coord(axis="x").coord_system = crs
        src_cube.coord(axis="y").coord_system = crs
        tgt_cube = ants.tests.stock.geodetic((4, 4))
        regridder = TwoStage()
        re = regridder.regridder(src_cube, tgt_cube)
        re(src_cube)

        self.assertTrue(self.area_patch.called)
        self.assertFalse(self.regtar_patch.called)
        self.assertFalse(self.linear_patch.called)


class Test___repr__(ants.tests.TestCase):
    def test(self):
        mdtol = 1
        scheme = TwoStage(mdtol=mdtol)
        tar = "TwoStage(mdtol={})".format(mdtol)
        self.assertEqual(repr(scheme), tar)


if __name__ == "__main__":
    ants.tests.main()
