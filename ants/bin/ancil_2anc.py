#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
"""
Conversion to ancillary file format application
***********************************************

Loads one or more cubes from the provided filepath(s) and saves them to the
specified output file path as both NetCDF and ancillary.  See
:ref:`input_output` for details on supported fileformats.

This is a file format translation utility. It assumes all fields in the input
file are to be added to the output file.  It will also infer most of the
metadata for the ancillary file format from the metadata in the input files.
To do this the metadata in the input files has to be complete and accurate.
This is discussed in more detail in :doc:`/appendixA_time_handling`.

As ancil_2anc is only a file format translation utility you may need to edit
the metadata of your input file to make it complete and accurate. You will
often need to specify the grid staggering as this usually can not be
inferred from the metadata in input files.

For an example usage of `ancil_2anc.py`, see :doc:`/examples`.

See Also
--------

`UM input and output file formats (F03) \
<https://code.metoffice.gov.uk/doc/um/latest/papers/umdp_F03.pdf>`_

"""
import warnings

import ants
import ants.io.save as save
import iris


def load_data(source):
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", "NetCDF default loading", iris._deprecation.IrisDeprecation
        )
        cubes = ants.load(source)
    return cubes


def main(source_path, output_path, grid_staggering, use_new_saver, netcdf_only):
    """
    Convert specified source file to an ancillary.

    Parameters
    ----------
    source_path : str
        Source data file path to be converted to an ancillary.
    output_path : str
        Output ancillary file path.  A NetCDF will also be written.
    grid_staggering : :obj:`int`, optional
        Grid staggering lookup, see UM F03 - FLH(9).

    Returns
    -------
    : :class:`~iris.cube.CubeList`
        The cubelist (which may contain one or more cubes) which has been
        written to an ancillary.

    Note
    ----
    Note that ANTS uses v2.3 of iris so CubeLists use ``extract`` with the
    ``strict`` argument rather than ``extract_cube``.

    """
    source_cubes = load_data(source_path)
    if grid_staggering is not None:
        for source_cube in source_cubes:
            source_cube.attributes["grid_staggering"] = grid_staggering

    if use_new_saver:
        if not netcdf_only:
            save.ancil(source_cubes, output_path)
        save.netcdf(source_cubes, output_path)
    else:
        ants.save(source_cubes, output_path)

    return source_cubes


def _get_parser():
    parser = ants.AntsArgParser()
    help_text = (
        "Specify the grid staggering of the source file, usually 3 "
        "(meaning New Dynamics) or 6 (meaning ENDGame).  This is "
        "needed if it is not present in the source (e.g.  PP "
        "sources) or is incorrect.  Translation between grid "
        "staggerings is not provided."
    )
    parser.add_argument("--grid-staggering", "-g", type=int, help=help_text)
    return parser


if __name__ == "__main__":
    args = _get_parser().parse_args()
    main(
        args.sources,
        args.output,
        args.grid_staggering,
        args.use_new_saver,
        args.netcdf_only,
    )
