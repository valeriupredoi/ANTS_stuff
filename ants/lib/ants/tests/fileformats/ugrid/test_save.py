# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import unittest.mock as mock
import warnings

import ants
import ants.tests
import iris
from ants.fileformats._ugrid import save


@mock.patch("ants.fileformats._ugrid.iris_save")
@mock.patch("ants.fileformats._ugrid._save_nodes")
@mock.patch("ants.fileformats._ugrid._save_minimal_topology")
@mock.patch("ants.fileformats._ugrid._save_UGRID_compliance")
class TestSaveCubeLists(ants.tests.TestCase):
    # Three cases: single cube, cubelist of a single cube, and a cubelist of
    # multiple cubes.  We accept all three.

    def test_single_cube(self, *args):
        cube = ants.tests.stock.mesh_C4()
        # Can't patch Dataset in class decorators since it's used in the load
        # in stock.mesh_C4():
        with mock.patch("ants.fileformats._ugrid.Dataset"):
            self.assertIsNone(save(cube, "foo"))

    def test_cubelist_of_single_cube(self, *args):
        cubes = iris.cube.CubeList([ants.tests.stock.mesh_C4()])
        # Can't patch Dataset in class decorators since it's used in the load
        # in stock.mesh_C4():
        with mock.patch("ants.fileformats._ugrid.Dataset"):
            self.assertIsNone(save(cubes, "foo"))

    def test_cubelist_of_multiple_cubes(self, *args):
        cubes = iris.cube.CubeList(
            [ants.tests.stock.mesh_C4(), ants.tests.stock.mesh_C4()]
        )
        cubes[1].rename("sample_data2")
        # Can't patch Dataset in class decorators since it's used in the load
        # in stock.mesh_C4():
        with mock.patch("ants.fileformats._ugrid.Dataset"):
            self.assertIsNone(save(cubes, "foo"))


@mock.patch("ants.fileformats._ugrid.iris_save")
@mock.patch("ants.fileformats._ugrid._save_nodes")
@mock.patch("ants.fileformats._ugrid._save_minimal_topology")
@mock.patch("ants.fileformats._ugrid._save_UGRID_compliance")
class TestExceptions(ants.tests.TestCase):
    def test_saver_argument_triggers_warning(self, *args):
        expected = warnings.warn(
            "Invalid fileformat foo specified.  Saving UGrid cubes "
            "bar as UGrid format."
        )

        cube = ants.tests.stock.mesh_C4()
        ants.config.CONFIG["saver"] = "foo"
        cube.rename("bar")

        # Can't patch Dataset in class decorators since it's used in the load.
        with mock.patch("ants.fileformats._ugrid.Dataset"):
            with self.assertWarnsRegex(UserWarning, expected):
                save(cube, "foo")


@mock.patch("ants.fileformats._ugrid.iris_save")
@mock.patch("ants.fileformats._ugrid.additional_save")
class TestXIOSWorkaround(ants.tests.TestCase):
    def test_workaround_called(self, *args):
        cube = ants.tests.stock.mesh_C4()
        # Can't patch Dataset in class decorators since it's used in the load
        # in stock.mesh_C4():
        with mock.patch(
            "ants.fileformats._ugrid._apply_XIOS_workaround"
        ) as mock_workaround:
            save(cube, "foo")
        mock_workaround.assert_called_once_with(iris.cube.CubeList([cube]))


if __name__ == "__main__":
    ants.tests.main()
