# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
from unittest import mock

import ants.tests
from ants.utils._dask import is_lazy_data


class TestAll(ants.tests.TestCase):
    def test_check_iris_is_lazy_data_spec(self):
        # Ensure private iris is_lazy_data function exists and matches the UI
        # we expect: We have to do this because it's private.
        patch = mock.patch("iris._lazy_data.is_lazy_data", spec_set=True)
        with patch as patched:
            is_lazy_data(mock.sentinel.data)
        self.assertTrue(patched.called)


if __name__ == "__main__":
    ants.tests.main()
