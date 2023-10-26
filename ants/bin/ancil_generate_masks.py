#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
"""
Generate masks application
**************************

Computes the land sea mask and sea mask from the land cover fraction, where
the land cover fraction typically comes from an ocean model.  Note that if
you're generating masks from the vegetation fractions, you should use the
`ancil_lct.py` application rather than this one.  For convenience, the land
cover fraction is also saved.  This will be equivalent to the source file.
Here are the steps taken:

- Extract the land cover fraction field ('m01s00i505').
- Derives the land sea mask (qrparm.mask) from the land cover type fraction
  field, where any point with land_frac > 0 is set as land (i.e. data value of
  1) in the qrparm.mask file.
- Derives the sea mask (qrparm.mask_sea) from the land cover type fraction
  field, where any point with land_frac < 1 is set as sea (i.e. data value of
  0) in the qrparm.mask file.

Fields returned:

- Land cover fraction, land mask and sea mask fields in that order
  ('m01s00i505', 'm01s00i030', 'm01s00i030' - note that both masks uses the
  same 'm01s00i030' code).

Note the ``output`` argument defines the directory where the mask files will
be written, rather than the paths for the individual masks.

"""

import os

import ants
import ants.decomposition
import ants.io.save as save
import iris
import numpy as np


def _prepare_mask_cube(cube):
    """
    Sets the attributes needed for mask files, and converts datatype to integer.

    Operates in place.

    Attributes 'valid_min' and 'valid_max' are set to ensure netCDF result is
    correct for boolean data.

    STASH code is set to m01s00i030.

    Data type is set to integer.

    Parameters
    ----------
    cube : :class:`iris.cube.Cube`
        Mask cube to prepare.

    Returns
    -------
    : None
    Operates in place.

    """
    cube.attributes["valid_min"] = 0
    cube.attributes["valid_max"] = 1
    cube.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi("m01s00i030")
    cube.data = cube.data.astype("int8")


def load_data(
    source_path,
):
    land_fraction_constraint = iris.AttributeConstraint(STASH="m01s00i505")
    source = ants.load_cube(source_path, constraint=land_fraction_constraint)
    return source


def _derive_masks(land_fraction):
    """
    Return land and sea masks derived from land fraction.

    For the land mask, any land fraction values > 0 are interpreted as land.
    For the sea mask, any land fraction value < 1 are interpreted as sea.

    Parameters
    ----------
    land_fraction : :class:`iris.cube.Cube`
        Land cover fraction cube

    Returns
    -------
    : tuple(:class:`iris.cube.Cube`, :class:`iris.cube.Cube`)
        The land mask and sea mask cubes, in that order.

    """
    # Coastal points (i.e. where land fraction is between 0 and 1) are classed
    # as land for the land mask and sea for the sea mask.  So we use a default
    # land mask value of 0 i.e. sea:
    land_mask = land_fraction.copy(data=np.zeros_like(land_fraction.data))
    # And where there's any land in the cell, we set the land mask to land:
    land_mask.data[land_fraction.data > 0] = 1

    # Similarly, default sea mask value is 1 i.e. land:
    sea_mask = land_fraction.copy(data=np.ones_like(land_fraction.data))
    # And where there's any sea in the cell, we set the sea mask to sea:
    sea_mask.data[land_fraction.data < 1] = 0

    land_mask.rename("land_binary_mask")
    sea_mask.rename("land_binary_mask")

    land_mask.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi("m01s00i030")
    sea_mask.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi("m01s00i030")

    return land_mask, sea_mask


def main(
    source_path,
    output_filepath,
    use_new_saver,
    netcdf_only,
):
    land_fraction = load_data(source_path)

    land_mask, sea_mask = ants.decomposition.decompose(_derive_masks, land_fraction)

    ants.config.dirpath_writeable(output_filepath)
    _prepare_mask_cube(land_mask)
    _prepare_mask_cube(sea_mask)

    cubes_to_save = {
        # "filename": (cube, fill_value),
        "qrparm.mask": (land_mask, -1),
        "qrparm.mask_sea": (sea_mask, -1),
        "qrparm.landfrac": (land_fraction, None),
    }

    for filename, (cube, fill_value) in cubes_to_save.items():
        filepath = os.path.join(output_filepath, filename)
        if use_new_saver:
            if not netcdf_only:
                save.ancil(cube, filepath)
            save.netcdf(cube, filepath, fill_value=fill_value)
        else:
            ants.save(cube, filepath, fill_value=fill_value)

    return land_fraction, land_mask, sea_mask


def _get_parser():
    parser = ants.AntsArgParser()
    return parser


if __name__ == "__main__":
    args = _get_parser().parse_args()
    main(
        args.sources,
        args.output,
        args.use_new_saver,
        args.netcdf_only,
    )
