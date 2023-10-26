# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import unittest.mock as mock

import ants.decomposition
import ants.tests
from ants.decomposition import DomainDecompose


class Test_src_generator(ants.tests.TestCase):
    def test_illdefined_relationship_1src_ntgt_alt_grid(self):
        # Ensure that an exception is raised if we have multiple targets, 1
        # source, yet the targets are not on idential grids as each other.
        decom = DomainDecompose()
        msg = "Ill-defined relationship between 1 source and multiple targets"
        target1 = ants.tests.stock.geodetic((1, 1))
        target2 = ants.tests.stock.geodetic((1, 2))
        mosaics = [
            ants.decomposition.MosaicBySplit(tgt, (1, 1)) for tgt in [target1, target2]
        ]
        decom._mosaics = mosaics
        decom._sources = [mock.sentinel.source]
        with self.assertRaisesRegex(RuntimeError, msg):
            decom.src_generator

    def test_illdefined_relationship_nsrc_ntgt_ambiguous_pairing(self):
        # Ambiguous pairing between multiple sources and multiple targets as
        # their number is greater than 1 but their number do not match.
        decom = DomainDecompose()
        msg = "Ill-defined relationship between number of sources and " "targets"
        decom._mosaics = [
            mock.sentinel.target1,
            mock.sentinel.target2,
            mock.sentinel.target3,
        ]
        decom._sources = [mock.sentinel.source1, mock.sentinel.source2]
        with self.assertRaisesRegex(RuntimeError, msg):
            decom.src_generator


if __name__ == "__main__":
    ants.tests.main()
