#!/usr/bin/env python2.7
"""Example to show filtering cubes by whether source data is a climatology
without relying on filename"""
import argparse
import glob

import iris


def load_cube_clim_noclim(filenames, clim):
    """
    Discriminate between climatology and non-climatology data by whether the
    time coordinate has 12 points or not.

    """
    # Because data is loaded lazily, this is a sufficiently efficient
    # solution.
    cubes = iris.load(glob.glob(in_files))
    if clim:
        cubes = [cube for cube in cubes if len(cube.coord(axis="t").points) == 12]
    else:
        cubes = [cube for cube in cubes if len(cube.coord(axis="t").points) != 12]

    cube = iris.cube.CubeList(cubes).merge_cube()
    return cube


#
# if __name__ == '__main__':
#    parser = argparse.ArgumentParser(
#        description=__doc__,
#        formatter_class=argparse.RawDescriptionHelpFormatter,
#    )
#    parser.add_argument(
#        '-c', '--clim', help='Load climatology data', action='store_true')
#    args = parser.parse_args()
#
#    in_files = ("/group_workspaces/jasmin2/tids/CMIP6_ANCIL/data/"
#                "inputs4MIPs_2017-10-10/UReading/surfaceFluxes/CMIP/"
#                "NCAR-CCMI-2-0/mon/drynoy/gn/v20161207/*")
#
#    cube = load_cube(in_files, args.clim)
#    print(cube.summary())
#
