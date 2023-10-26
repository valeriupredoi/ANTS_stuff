#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
"""
General regrid application
**************************

Regrids data from a source to a target grid using
:class:`ants.regrid.GeneralRegridScheme`.  The result is written to an output
file.  The application supports both horizontal and vertical regridding.  The
regrid algorithm can be specified in the ants configuration file as described
in :class:`ants.config.GlobalConfiguration`. See :mod:`ants.regrid` and
:doc:`/regrid` for further
details.

If a landseamask is provided using the `target-lsm` argument, this mask is
honoured by only filling missing land points with valid land values (or
similarly, missing sea points with valid sea values).

For example usages of `ancil_general_regrid.py`, see :doc:`/examples`.

"""
import warnings

import ants
import ants.decomposition as decomp
import ants.io.save as save
import ants.utils
from ants.utils.cube import create_time_constrained_cubes


def load_data(
    source,
    target_grid=None,
    target_landseamask=None,
    land_fraction_threshold=None,
    begin=None,
    end=None,
):
    source_cubes = ants.load(source)
    if begin is not None:
        source_cubes = create_time_constrained_cubes(source_cubes, begin, end)
    if target_grid:
        target_cube = ants.load_grid(target_grid)
    else:
        target_cube = ants.fileformats.load_landsea_mask(
            target_landseamask, land_fraction_threshold
        )
    return source_cubes, target_cube


def regrid(sources, target):
    sources = ants.utils.cube.as_cubelist(sources)
    results = []
    scheme = ants.regrid.GeneralRegridScheme()
    if scheme._horizontal_scheme is None:
        scheme = ants.regrid.GeneralRegridScheme(horizontal_scheme="TwoStage")
        msg = (
            "No horizontal regridding scheme has been set and so the "
            "TwoStage regridder is being used. The default setting of "
            "the horizontal regridder scheme to TwoStage has been "
            "deprecated as of ANTS 1.1.0 and will be removed in a "
            "future release. Please set a horizontal regridding scheme in "
            "the ANTS config."
        )
        warnings.warn(msg, FutureWarning, stacklevel=2)
    for source in sources:
        results.append(source.regrid(target, scheme))
    return results


def main(
    source_path,
    output_path,
    target_path,
    target_lsm_path,
    invert_mask,
    land_fraction_threshold,
    begin,
    end,
    save_ukca,
    use_new_saver,
    netcdf_only,
):
    """
    General regrid application top level call function.

    Loads source data cubes, regrids them to match target data cube
    co-ordinates, and saves result to output.  In addition to writing the
    resulting data cube to disk, also returns the regridded data cube.

    Parameters
    ----------

    source_path : str
        File path for one or more files which contain the data to be
        regridded.
    target_path : str
        File path for files that provide the grid to which the source data
        cubes will be mapped.  Separate files can be provided to generate a
        complete grid i.e. a namelist for vertical levels can be used with a
        data file for the horizontal coordinates.
    target_lsm_path : str
        File path for a land sea mask that provides the grid to which
        the source data cube will be mapped.  The output will be made
        consistent with this land sea mask.
    invert_mask : :obj:`bool`, optional
        Determines whether to invert the mask for the `target_lsm_path`
        argument.
        When set to True, treat target mask True (1) values as unmasked.  When set
        to False, treat target mask True (1) values as masked. Default is True.
    output_path : str
        Output file path to write the regridded data to.
    land_fraction_threshold : str
    begin : :obj:`datetime`, optional
        If provided, all source data prior to this year is discarded.  Default is to
        include all source data.
    end : :obj:`datetime`, optional
        If provided, all source data after this year is discarded.  Default is to
        include all source data.
    Returns
    -------
    : :class:`~iris.cube.Cube`
    A single data cube with the regridded data.

    """
    source_cubes, target_cube = load_data(
        source_path,
        target_path,
        target_lsm_path,
        land_fraction_threshold,
        begin,
        end,
    )
    if ants.utils.cube._is_ugrid(target_cube):
        raise ValueError(
            "Target appears to be a UGrid mesh - the regrid to mesh application in "
            "ANTS contrib should be used instead."
        )

    regridded_cubes = decomp.decompose(regrid, source_cubes, target_cube)
    if target_lsm_path:
        ants.analysis.make_consistent_with_lsm(
            regridded_cubes, target_cube, invert_mask
        )

    if use_new_saver:
        if save_ukca:
            save.ukca_netcdf(regridded_cubes, output_path)
        else:
            if not netcdf_only:
                save.ancil(regridded_cubes, output_path)
            save.netcdf(regridded_cubes, output_path)
    else:
        ants.save(regridded_cubes, output_path)

    return regridded_cubes


def _get_parser():
    parser = ants.AntsArgParser(
        target_lsm=True, target_grid=True, time_constraints=True
    )
    invmask_help = (
        "Invert the provided target_mask or not.\n"
        "Using this argument will set it to False. "
        "When set to True, treat target mask True (1) values as unmasked. When set "
        "to False, treat target mask True (1) values as masked. "
        "It is common to use this argument to denote source ocean fields as the "
        "landsea mask has True values to denote land."
    )
    parser.add_argument(
        "--invert-mask",
        action="store_false",
        help=invmask_help,
        required=False,
    )
    parser.add_argument(
        "--save-ukca",
        action="store_true",
        help="Save to a UKCA-specific netCDF file",
        required=False,
    )
    return parser


if __name__ == "__main__":
    parser = _get_parser()
    args = parser.parse_args()

    source = args.sources
    main(
        source,
        args.output,
        args.target_grid,
        args.target_lsm,
        args.invert_mask,
        args.land_threshold,
        args.begin,
        args.end,
        args.save_ukca,
        args.use_new_saver,
        args.netcdf_only,
    )
