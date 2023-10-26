# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
from ants.utils.cube import is_equal_hgrid


class TestRegular(ants.tests.TestCase):
    def test_equal(self):
        cube1 = ants.tests.stock.geodetic((1, 1))
        cube2 = ants.tests.stock.geodetic((1, 1))
        self.assertTrue(is_equal_hgrid([cube1, cube2]))

    def test_not_equal(self):
        cube1 = ants.tests.stock.geodetic((1, 1))
        cube2 = ants.tests.stock.geodetic((1, 1), ylim=(-80, 80))
        self.assertFalse(is_equal_hgrid([cube1, cube2]))


class TestUGrid(ants.tests.TestCase):
    def test_mixed_cubes(self):
        cube1 = ants.tests.stock.geodetic((1, 1))
        cube2 = ants.tests.stock.mesh_C4()

        self.assertFalse(is_equal_hgrid([cube1, cube2]))

    def test_different_coords(self):
        cube1 = ants.tests.stock.mesh_C4()
        cube2 = ants.tests.stock.mesh_C4()
        cube2.coord("latitude").points += 0.1

        self.assertFalse(is_equal_hgrid([cube1, cube2]))

    def test_node_attributes(self):
        cube1 = ants.tests.stock.mesh_C4()
        cube2 = ants.tests.stock.mesh_C4()
        cube2.attributes["nodes"] = "foo"

        self.assertFalse(is_equal_hgrid([cube1, cube2]))

    def test_different_mapping(self):
        cube1 = ants.tests.stock.mesh_C4()
        cube2 = ants.tests.stock.mesh_C4()
        cube2.attributes["face_node_connectivity"] = "foo"

        self.assertFalse(is_equal_hgrid([cube1, cube2]))

    def test_equal(self):
        cube1 = ants.tests.stock.mesh_C4()
        cube2 = ants.tests.stock.mesh_C4()

        self.assertTrue(is_equal_hgrid([cube1, cube2]))

    def test_absent_attribute_on_all_cubes_is_equal(self):
        cube1 = ants.tests.stock.mesh_C4()
        cube2 = ants.tests.stock.mesh_C4()
        [cube.attributes.pop("face_node_connectivity") for cube in (cube1, cube2)]

        self.assertTrue(is_equal_hgrid([cube1, cube2]))

    def test_absent_attribute_on_some_cubes_is_not_equal(self):
        cube1 = ants.tests.stock.mesh_C4()
        cube2 = ants.tests.stock.mesh_C4()

        cube2.attributes.pop("face_node_connectivity")

        self.assertFalse(is_equal_hgrid([cube1, cube2]))


if __name__ == "__main__":
    ants.tests.main()
