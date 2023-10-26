# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import unittest

import ants.tests
import dask.array as da
from ants.utils._dask import copy


class DummyDataProxy(object):
    def __init__(self, shape, dtype):
        self.shape = shape
        self.dtype = dtype
        self.recorded_keys = []

    @property
    def ndim(self):
        return len(self.shape)

    def __getitem__(self, keys):
        self.recorded_keys.append(keys)
        return None


class Test_concatenation(ants.tests.TestCase):
    # Dask gives undesrable behaviour (chunking) when it comes to
    # concatenating dask arrays, see
    # https://code.metoffice.gov.uk/trac/ancil/ticket/1098
    def setUp(self):
        self.proxy = DummyDataProxy((32, 20, 30), "float64")
        self.darray = da.from_array(self.proxy, chunks=(16, 20, 30))

    def check_concat_chunking(self, da1, da2):
        res = da.concatenate([da1, da2], 2)
        try:
            res.compute()
        except IndexError:
            pass
        tar = [
            (slice(0, 0, None), slice(0, 0, None), slice(0, 0, None)),
            (slice(0, 0, None), slice(0, 0, None), slice(0, 0, None)),
            (slice(0, 0, None), slice(0, 0, None), slice(0, 0, None)),
            (slice(0, 16, None), slice(0, 4, None), slice(28, 30, None)),
            (slice(16, 32, None), slice(0, 4, None), slice(28, 30, None)),
            (slice(0, 16, None), slice(0, 4, None), slice(0, 2, None)),
            (slice(16, 32, None), slice(0, 4, None), slice(0, 2, None)),
        ]
        self.assertTrue(len(self.proxy.recorded_keys), 7)
        self.assertEqual(sorted(self.proxy.recorded_keys), sorted(tar))

    @unittest.expectedFailure
    def test_native(self):
        # If this should suceed, we would then likely no longer need the
        # safeslicing workaround anymore.
        da1 = self.darray[:, :4, :2]
        da2 = self.darray[:, :4, -2:]
        self.check_concat_chunking(da1, da2)

    def test_workaround(self):
        da1 = copy(self.darray)[:, :4, :2]
        da2 = copy(self.darray)[:, :4, -2:]
        self.check_concat_chunking(da1, da2)


if __name__ == "__main__":
    ants.tests.main()
