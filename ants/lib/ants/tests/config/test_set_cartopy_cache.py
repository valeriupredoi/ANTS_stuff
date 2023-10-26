# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import os
from unittest import mock

import ants.tests
import cartopy
from ants.config import set_cartopy_cache


@mock.patch.dict(os.environ, {"ANTS_CARTOPY_CACHE": "my/fake/path"})
class TestSetCartopyCacheLocationByEnvVar(ants.tests.TestCase):
    # cartopy is not re-imported for each test. Hence, need to restore the
    # original default value for the data_dir key after this test, to
    # prevent future tests creating a `my/fake/path` directory
    def test_set_cache_location(self):
        original_data_dir = cartopy.config["data_dir"]
        set_cartopy_cache()
        self.assertEqual(cartopy.config["data_dir"], "my/fake/path")
        cartopy.config["data_dir"] = original_data_dir

    def test_logging_of_cache_location(self):
        original_data_dir = cartopy.config["data_dir"]
        with self.assertLogs("ants.config", level="INFO") as context_manager:
            set_cartopy_cache()
            self.assertEqual(
                context_manager.output,
                ["INFO:ants.config:ANTS_CARTOPY_CACHE is set to: my/fake/path"],
            )
        cartopy.config["data_dir"] = original_data_dir


class TestCartopyCacheLocationSetByConfig(ants.tests.TestCase):
    def test_not_setting_cache_location_causes_no_change(self):
        expected = cartopy.config["data_dir"]

        with mock.patch.dict(os.environ):
            set_cartopy_cache()
            self.assertEqual(cartopy.config["data_dir"], expected)


if __name__ == "__main__":
    ants.tests.main()
