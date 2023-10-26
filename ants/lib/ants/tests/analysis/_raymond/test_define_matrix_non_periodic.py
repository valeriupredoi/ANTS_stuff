# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
import numpy as np
from ants.analysis._raymond import define_matrix_non_periodic as define_matrix


class TestAll(ants.tests.TestCase):
    # Checking the return values of the matrix initialisation.
    def test_values(self):
        res = define_matrix(7, 1)
        tar = np.array(
            [
                [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [2.0, 0.0, 12.0, 0.0, 2.0, 0.0, 0.0],
                [0.0, 12.0, 0.0, 40.0, 0.0, 12.0, 0.0],
                [0.0, 0.0, 2.0, 0.0, 12.0, 0.0, 2.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
            ]
        )
        self.assertArrayEqual(res.toarray(), tar)


if __name__ == "__main__":
    ants.tests.main()
