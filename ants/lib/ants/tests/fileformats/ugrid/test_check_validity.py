# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
import iris
from ants.fileformats._ugrid import check_validity


class TestAll(ants.tests.TestCase):
    def test_exception_for_different_mesh_attributes(self):
        cube1 = ants.tests.stock.mesh_C4()
        cube2 = cube1.copy()
        cube1.rename("cube1")
        cube2.rename("cube2")
        cube2.attributes["face_node_connectivity"] = "foo"
        cubes = iris.cube.CubeList((cube1, cube2))

        with self.assertRaisesRegex(IOError, "Cannot save cubes with differ"):
            check_validity(cubes, "latitude", "longitude")


if __name__ == "__main__":
    ants.tests.main()
