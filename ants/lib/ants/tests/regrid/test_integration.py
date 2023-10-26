# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import unittest.mock as mock

import ants.tests
from ants.regrid import GeneralRegridScheme


class TestErrorMessages(ants.tests.TestCase):
    def test_no_scheme_given(self):
        source = ants.tests.stock.geodetic((2, 2))
        target = ants.tests.stock.geodetic((2, 2))
        scheme = GeneralRegridScheme()
        with self.assertRaises(AttributeError) as context:
            source.regrid(target, scheme)

            self.assertTrue(
                "At least one of horizontal \
            or vertical re-grid schemes must be provided."
                in context.exception
            )


class TestInterpolation(ants.tests.TestCase):
    def setUp(self):
        self.source = ants.tests.stock.geodetic((2, 2))
        self.target = self.source.copy()

    def test_conservative(self):
        scheme = "ants.regrid.interpolation.Conservative"
        with mock.patch(scheme) as patched_scheme:
            scheme = GeneralRegridScheme(vertical_scheme="Conservative")
            self.source.regrid(self.target, scheme)
        self.assertTrue(patched_scheme.called)

    def test_linear(self):
        scheme = "ants.regrid.interpolation.Linear"
        with mock.patch(scheme) as patched_scheme:
            scheme = GeneralRegridScheme(vertical_scheme="Linear")
            self.source.regrid(self.target, scheme)
        self.assertTrue(patched_scheme.called)

    def test_nearest(self):
        scheme = "ants.regrid.interpolation.Nearest"
        with mock.patch(scheme) as patched_scheme:
            scheme = GeneralRegridScheme(vertical_scheme="Nearest")
            self.source.regrid(self.target, scheme)
        self.assertTrue(patched_scheme.called)


class TestESMF(ants.tests.TestCase):
    @ants.tests.skip_esmf
    def test_conservative(self):
        scheme = GeneralRegridScheme(horizontal_scheme="ConservativeESMF")
        source = ants.tests.stock.geodetic((2, 2))
        target = source.copy()
        res = source.regrid(target, scheme)
        self.assertEqual(res, target)

    def test_expected_scheme(self):
        with mock.patch("ants.regrid.esmf.ConservativeESMF") as patched_scheme:
            scheme = GeneralRegridScheme(horizontal_scheme="ConservativeESMF")
            source = ants.tests.stock.geodetic((2, 2))
            target = source.copy()
            source.regrid(target, scheme)
        self.assertTrue(patched_scheme.called)


class TestRectilinear(ants.tests.TestCase):
    def test_twostage(self):
        with mock.patch("ants.regrid.rectilinear.TwoStage") as patched_scheme:
            scheme = GeneralRegridScheme(horizontal_scheme="TwoStage")
            source = ants.tests.stock.geodetic((2, 2))
            target = source.copy()
            source.regrid(target, scheme)
        self.assertTrue(patched_scheme.called)

    def test_expected_scheme(self):
        # Ensure our Linear is used over the iris one.
        with mock.patch("ants.regrid.rectilinear.Linear") as patched_scheme:
            scheme = GeneralRegridScheme(horizontal_scheme="Linear")
            source = ants.tests.stock.geodetic((2, 2))
            target = source.copy()
            source.regrid(target, scheme)
        self.assertTrue(patched_scheme.called)


class TestIris(ants.tests.TestCase):
    def test_areaweighted(self):
        scheme = GeneralRegridScheme(horizontal_scheme="AreaWeighted")
        source = ants.tests.stock.geodetic((2, 2))
        target = source.copy()
        res = source.regrid(target, scheme)
        self.assertEqual(res, target)

    def test_expected_scheme(self):
        with mock.patch("ants.regrid.rectilinear.AreaWeighted") as patched_scheme:
            scheme = GeneralRegridScheme(horizontal_scheme="AreaWeighted")
            source = ants.tests.stock.geodetic((2, 2))
            target = source.copy()
            source.regrid(target, scheme)
        self.assertTrue(patched_scheme.called)


if __name__ == "__main__":
    ants.tests.main()
