# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
from ants.fileformats import set_saver
from iris.io import save


class TestDetermineFormat(ants.tests.TestCase):
    def test_set_saver(self):
        savers = ["nc", "pp"]
        for saver in savers:
            with self.subTest():
                target = "test_filename"
                result = set_saver(saver, target)
                self.assertTrue(callable(result))
                self.assertEqual(result.__name__, "save")

    def test_set_ancil_saver(self):
        target = "test_filename"
        result = set_saver("ancil", target)
        self.assertTrue(callable(result))
        self.assertEqual(result.__name__, "_enforced_netcdf_on_ancil_save")

    def test_unrecognised_file_raises_error(self):
        saver_spec = "blah"
        target = "test_filename"
        with self.assertRaises(ValueError):
            set_saver(saver_spec, target)
            self.assertRaisesRegex(
                ValueError,
                'Cannot save; no saver can be found associated with "blah"',
                save,
            )


if __name__ == "__main__":
    ants.tests.main()
