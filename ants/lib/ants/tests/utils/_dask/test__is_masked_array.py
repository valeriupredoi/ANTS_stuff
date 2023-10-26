# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
import dask
import numpy as np
from ants.utils._dask import _is_masked_array


class TestAll(ants.tests.TestCase):
    def test__is_masked_array_for_masked_array(self):
        masked_array = np.ma.array([1, 2, 3], mask=[0, 1, 0])
        dask_array = dask.array.from_array(masked_array)
        # Actually checking the dask behaviour here as we are using a private
        # attribute that is undocumented (_meta).
        self.assertTrue(_is_masked_array(dask_array))

    def test__is_masked_array_for_unmasked_array(self):
        unmasked_array = np.array([1, 2, 3])
        dask_array = dask.array.from_array(unmasked_array)
        # Actually checking the dask behaviour here as we are using a private
        # attribute that is undocumented (_meta).
        self.assertFalse(_is_masked_array(dask_array))


if __name__ == "__main__":
    ants.tests.main()
