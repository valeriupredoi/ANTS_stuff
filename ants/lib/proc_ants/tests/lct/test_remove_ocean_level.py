# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
import iris
import numpy as np
from proc_ants.lct import remove_ocean_level


class TestAll(ants.tests.TestCase):
    def test_all(self):
        data = np.array([[[0.2], [0.6]], [[0.4], [0.2]], [[0.4], [0.2]]])
        lct = ants.tests.stock.geodetic(data.shape, data=data)
        pseudo_coord = iris.coords.AuxCoord([0, 1, 2], long_name="pseudo_level")

        lct.add_aux_coord(pseudo_coord, 0)

        res, lsm = remove_ocean_level(lct)
        tar = np.ma.array([[[0.5], [0.5]], [[0.5], [0.5]]])
        tar[-1] = np.ma.masked
        lsm_tar = np.ma.array([[1], [0]])
        self.assertArrayAlmostEqual(res.data, tar)
        self.assertArrayEqual(lsm.data, lsm_tar)


if __name__ == "__main__":
    ants.tests.main()
