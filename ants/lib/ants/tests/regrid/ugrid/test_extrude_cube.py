# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
import iris
import numpy as np
from ants.regrid._ugrid import _extrude_cube as extrude_cube


class TestAll(ants.tests.TestCase):
    def setUp(self):
        self.cube = ants.tests.stock.mesh_C4()

    def test_returns_cube(self):
        expected = iris.cube.Cube

        actual = extrude_cube(self.cube, (3, 96), dtype=self.cube.dtype)

        self.assertIsInstance(actual, expected)

    def test_result_shape(self):
        expected = (3, 96)

        actual = extrude_cube(self.cube, expected, dtype=self.cube.dtype).shape

        self.assertEqual(expected, actual)

    def test_reject_invalid_shape(self):
        with self.assertRaisesRegex(ValueError, "Invalid shape"):
            extrude_cube(self.cube, (3, 2), dtype=self.cube.dtype)

    def test_data_initialised_to_zero(self):
        expected = np.zeros((3, 96))

        actual = extrude_cube(self.cube, (3, 96), dtype=self.cube.dtype).data

        self.assertArrayEqual(expected, actual)

    def test_data_type(self):
        expected = self.cube.dtype

        actual = extrude_cube(self.cube, (3, 96), dtype=self.cube.dtype).dtype

        self.assertEqual(expected, actual)

    def test_long_name(self):
        self.cube.long_name = "foo"
        expected = self.cube.long_name

        actual = extrude_cube(self.cube, (3, 96), dtype=self.cube.dtype).long_name

        self.assertEqual(expected, actual)

    def test_attributes(self):
        self.cube.attributes["foo"] = "bar"
        expected = self.cube.attributes

        actual = extrude_cube(self.cube, (3, 96), dtype=self.cube.dtype).attributes

        self.assertEqual(expected, actual)

    def test_units(self):
        self.cube.units = "m"
        expected = self.cube.units

        actual = extrude_cube(self.cube, (3, 96), dtype=self.cube.dtype).units

        self.assertEqual(expected, actual)

    def test_horizontal_coordinates(self):
        expected = self.cube.coords(dimensions=0)

        result = extrude_cube(self.cube, (3, 96), dtype=self.cube.dtype)
        actual = result.coords(dimensions=1)

        self.assertEqual(expected, actual)

    def test_original_cube_unchanged(self):
        expected = self.cube.copy()

        extrude_cube(self.cube, (3, 96), dtype=self.cube.dtype)
        actual = self.cube

        self.assertEqual(expected, actual)


if __name__ == "__main__":
    ants.tests.main()
