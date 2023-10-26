#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
"""
Add lakes application
*********************

Adds the OSTIA lakes from preprocessed CCI data to the land sea mask.  Here
are the steps taken:

- Extracts the flag for resolved lake from the preprocessed CCI data.
- Performs a nearest neighbour regrid to get the CCI data onto the same grid
  as the land sea mask.
- Sets mask points as ocean where the corresponding vegetation fraction point
  is resolved lake.

Fields returned:

- Land sea mask ('m01s00i030')

"""
import ants
import ants.io.save as save
import iris


def load_data(source_path, vegetation_fraction_path):
    source = ants.load_cube(source_path)
    vegetation_fraction = ants.load_cube(vegetation_fraction_path)
    return source, vegetation_fraction


def _get_resolved_lake_value(vegetation_fraction):
    flag_meanings = vegetation_fraction.attributes["flag_meanings"].split()
    flag_values = vegetation_fraction.attributes["flag_values"]
    resolved_lake_index = flag_meanings.index("resolved_lake")
    resolved_lake_value = flag_values[resolved_lake_index]

    return resolved_lake_value


def _add_lakes(source, vegetation_fraction, resolved_lake_value):
    """
    Return land and sea masks where OSTIA lakes have been added.

    Parameters
    ----------
    source :  :class:`iris.cube.Cube`
        Land sea mask cube
    vegetation_fraction : :class:`iris.cube.Cube`
        Land cover fraction cube

    Returns
    -------
    : :class:`iris.cube.Cube`
        Source cube with OSTIA lakes added from the vegetation_fraction

    """

    regridded_vegetation_fraction = vegetation_fraction.regrid(
        source, scheme=iris.analysis.Nearest()
    )
    source.data[regridded_vegetation_fraction.data == resolved_lake_value] = 0

    return source


def main(
    source_path,
    vegetation_path,
    output_filepath,
    use_new_saver,
    netcdf_only,
):
    frac, vegetation_fraction = load_data(source_path, vegetation_path)

    resolved_lake_value = _get_resolved_lake_value(vegetation_fraction)
    result = _add_lakes(frac, vegetation_fraction, resolved_lake_value)

    ants.config.dirpath_writeable(output_filepath)

    if use_new_saver:
        if not netcdf_only:
            save.ancil(frac, output_filepath)
        save.netcdf(frac, output_filepath, fill_value=-1)
    else:
        ants.save(frac, output_filepath, fill_value=-1)

    return result


def _get_parser():
    parser = ants.AntsArgParser()
    parser.add_argument(
        "--vegetation-fraction",
        type=str,
        help="Vegetation fraction file defining lakes to use",
        required=True,
    )
    return parser


if __name__ == "__main__":
    args = _get_parser().parse_args()
    main(
        args.sources,
        args.vegetation_fraction,
        args.output,
        args.use_new_saver,
        args.netcdf_only,
    )
