# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
import iris
import numpy as np
from proc_ants.lct import set_whole_fraction_ice


class TestAll(ants.tests.TestCase):
    def test_all(self):
        data = np.array(
            [[[0.1], [0.2]], [[0.1], [0.2]], [[0.1], [0.4]], [[0.7], [0.2]]]
        )
        lct = ants.tests.stock.geodetic(data.shape, data=data)
        pseudo_coord = iris.coords.AuxCoord([0, 1, 2, 9], long_name="pseudo_level")

        lct.add_aux_coord(pseudo_coord, 0)

        set_whole_fraction_ice(lct)
        tar = [[[0.0], [0.25]], [[0.0], [0.25]], [[0.0], [0.5]], [[1.0], [0.0]]]
        self.assertArrayAlmostEqual(lct.data, tar)


if __name__ == "__main__":
    ants.tests.main()
