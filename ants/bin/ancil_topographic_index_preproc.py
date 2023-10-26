#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
"""
Pre-processing topographic index application
********************************************

Apply selected modifications to the base-line topographic index dataset, which
is supplied in netCDF format.

The following modifications are applied:

* For the topidx variable, set the fill value to -1.
* For the topidx variable, set the long_name attribute to 'topographic index'.
* Update the history attribute to reflect these changes.

"""
import ants
import ants.config
import ants.io.save as save
import ants.utils
import iris
import numpy as np


def load_data(filename):
    """Load the source data into an Iris cube."""
    constraint = iris.Constraint(cube_func=lambda cube: cube.var_name == "topidx")
    cube = ants.load_cube(filename, constraint)
    return cube


def process_data(cube):
    """Apply any pre-processing to the source data."""
    # Mask out data values equal to -1 by using that as the new fill value.
    if np.ma.isMA(cube.data):
        cube.data.set_fill_value(-1.0)
        cube.data[cube.data == -1] = np.ma.masked
    else:
        cube.data = np.ma.masked_values(cube.data, -1.0, copy=False)

    # Data is missing below a certain latitude.  This is remedied in the
    # topographic index generation application.
    ants.utils.cube.guess_horizontal_bounds(cube)

    bounds = cube.coord("latitude").bounds.copy()
    bounds[0][0] = -90
    bounds[-1][-1] = 90
    cube.coord("latitude").bounds = bounds
    bounds = cube.coord("longitude").bounds.copy()
    bounds[-1][-1] = 180
    cube.coord("longitude").bounds = bounds
    return cube


def main(input_filepath, output_filepath, use_new_saver):
    cube = load_data(input_filepath)
    cube = process_data(cube)
    cube.long_name = "topographic index"

    if use_new_saver:
        save.netcdf(cube, output_filepath)
    else:
        ants.save(cube, output_filepath)

    return cube


def _get_parser():
    parser = ants.AntsArgParser()
    return parser


if __name__ == "__main__":
    args = _get_parser().parse_args()
    main(args.sources, args.output, args.use_new_saver)
