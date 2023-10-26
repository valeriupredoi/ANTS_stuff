#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
"""
Pre-processor for the ITE (CEH) source for land cover types
***********************************************************

The following steps are taken:

 1. Set NetCDF "source" attribute metadata.
 2. Set the flag meanings and flag values for the dataset.

"""
import ants
import ants.config
import ants.io.save as save
import ants.utils
import proc_ants.lct_preproc_ite as lct_preproc_ite


def load_data(sources_filepath):
    cube = ants.utils.cube.concatenate_cube(ants.load(sources_filepath))
    return cube


def main(sources_filepath, output_filepath, use_new_saver):
    ite_cube = load_data(sources_filepath)
    lct_preproc_ite.pre_process_source(ite_cube)

    if use_new_saver:
        save.netcdf(ite_cube, output_filepath, fill_value=0)
    else:
        ants.save(ite_cube, output_filepath, fill_value=0)

    return ite_cube


def _get_parser():
    parser = ants.AntsArgParser()
    return parser


if __name__ == "__main__":
    args = _get_parser().parse_args()
    main(args.sources, args.output, args.use_new_saver)
