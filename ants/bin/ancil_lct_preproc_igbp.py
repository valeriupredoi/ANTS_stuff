#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
"""
Pre-processor for the IGBP source for land cover types
******************************************************

The following steps are taken:

 1. Introduce distinction between ocean and in-land water using the BATS source
    where flag value 15 is ocean in the BATS source (identical grid as the
    IGBP).
 2. Set missing/correct metadata:

    - Populate flag_values.
    - Populate flag_meanings.
    - Correct x and y bounds.
    - Populate missing coordinate system (WGS84).
    - Provide a suitable long_name (IGBP land classification).

 3. Fill missing data using the spiral search.

.. warning::

    The IGBP dataset is fraught with problems such as non-linear geolocation
    misalignment and missing islands.  It is recommended to use the global CCI
    (ESA) dataset instead as it is superior in every regard (resolution,
    accuracy, number of classification types etc.).

"""
import ants
import ants.config
import ants.io.save as save
import ants.utils
from proc_ants import lct_preproc_igbp


def load_data(igbp_filepath, bats_filepath):
    igbp_cube = ants.load_cube(igbp_filepath)
    bats_cube = ants.load_cube(bats_filepath)
    return igbp_cube, bats_cube


def main(igbp_filepath, bats_filepath, output_filepath, use_new_saver):
    igbp_cube, bats_cube = load_data(igbp_filepath, bats_filepath)
    lct_preproc_igbp.pre_process(igbp_cube, bats_cube)

    filler = ants.analysis.FillMissingPoints(igbp_cube)
    filler(igbp_cube)

    if use_new_saver:
        save.netcdf(igbp_cube, output_filepath, fill_value=100)
    else:
        ants.save(igbp_cube, output_filepath, fill_value=100)

    return igbp_cube


def _get_parser():
    parser = ants.AntsArgParser()
    parser.add_argument(
        "--bats-source",
        type=str,
        help=(
            "Pathname of BATS source, used to distinguish "
            "between igbp inland water and ocean"
        ),
        required=True,
    )
    return parser


if __name__ == "__main__":
    args = _get_parser().parse_args()
    main(args.sources, args.bats_source, args.output, args.use_new_saver)
