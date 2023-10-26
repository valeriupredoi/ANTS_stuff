# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
import iris
import numpy as np
from ants.regrid._ugrid import _build_mesh as build_mesh


@ants.tests.skip_esmf
class TestAll(ants.tests.TestCase):
    def test_type(self):
        import ESMF

        expected = ESMF.Mesh

        actual = build_mesh(ants.tests.stock.mesh_C4())

        self.assertIsInstance(actual, expected)

    def test_nodes_not_masked(self):
        # See
        # http://www.earthsystemmodeling.org/esmf_releases/public/ESMF_7_1_0r/esmpy_doc/html/mesh.html#ESMF.api.mesh.Mesh
        # for details of mesh.mask hardcoded indices
        source = ants.tests.stock.mesh_C4()
        source.data = np.ma.masked_array(source.data)
        source.data[0] = np.ma.masked
        mesh = build_mesh(source)
        actual = mesh.mask[0]

        self.assertIsNone(actual)

    def test_element_mask_2d(self):
        expected = np.zeros((96,), dtype=np.int32)
        expected[10] = 1

        source = ants.tests.stock.mesh_C4()
        source.data = np.ma.masked_array(source.data)
        source.data[10] = np.ma.masked
        mesh = build_mesh(source)
        actual = mesh.mask[1]

        self.assertArrayEqual(expected, actual)

    def test_element_mask_Nd(self):
        expected = np.zeros((96,), dtype=np.int32)
        expected[10] = 1

        source = ants.tests.stock.mesh_C4()
        source.data = np.ma.masked_array(np.arange(96))
        source.data[10] = np.ma.masked
        source2 = source.copy()
        source.add_aux_coord(iris.coords.AuxCoord(0, long_name="bing"), None)
        source2.add_aux_coord(iris.coords.AuxCoord(2, long_name="bing"), None)
        source2.data += 20
        source = iris.cube.CubeList([source, source2]).merge_cube()
        mesh = build_mesh(source)
        actual = mesh.mask[1]

        self.assertArrayEqual(expected, actual)


if __name__ == "__main__":
    ants.tests.main()
