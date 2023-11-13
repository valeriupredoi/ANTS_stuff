#!/usr/bin/env python
"""
Create the ancillaries needed to represent nitrogen deposition in UKESM1 for the CMIP6/DEC historical simulation 
00447    0 NITROGEN DEPOSITION (kgN/m2/s)
"""
import argparse
import datetime
import glob

import cf_units
import iris
import iris.analysis
import iris.fileformats
import iris.plot as iplt
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np
import pylab
from iris.time import PartialDateTime

import ants
import ants.decomposition as decomp
import ants.fileformats
import ants.io.save as save

print("Here")
iris.FUTURE.netcdf_promote = True

a4port = (8.27, 11.69)  # A4 paper portrait in inches
a4land = (11.69, 8.27)  # A4 paper landscape in inches
colour = ["green", "blue", "orange", "red", "purple", "black"]
#
# 30th January 2019: Revision 7236 of this vn26 script works, to the extent that I have tested it below.
#
# THese two command created historical ancillaries (a few years in lenght) whcih look OK and give same data as the (11/1895 looked at):
# /data/d04/hadsl/ANCILS/UKESM1/NITROGEN_DEPOSITION/n96e_rev5199/ndep_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-2-0_gn_184901-201512.n96e.anc.nc
#
# python3 cmip6_ndep_jasmin_vn26_past_and_future.py -s './TESTDIR_N96E/drynhx_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-2-0_gn_185001-201412.nc'
#   './TESTDIR_N96E/drynoy_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-2-0_gn_185001-201412.nc'
#   './TESTDIR_N96E/wetnhx_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-2-0_gn_185001-201412.nc'
#   './TESTDIR_N96E/wetnoy_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-2-0_gn_185001-201412.nc'
# --resolution  './TESTDIR_N96E_SSP585/qrparm.mask_n96e_eorca1_v2.2x_ESMF'
# -o './TESTDIR_JAN2019/Ndep_n96e_1893_1895_appended' -n 'Ndep' --begin 1893 --end 1895 -a true
#
# python3 cmip6_ndep_jasmin_vn26_past_and_future.py -s './TESTDIR_N96E/drynhx_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-2-0_gn_185001-201412.nc'
#   './TESTDIR_N96E/drynoy_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-2-0_gn_185001-201412.nc'
#   './TESTDIR_N96E/wetnhx_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-2-0_gn_185001-201412.nc'
#   './TESTDIR_N96E/wetnoy_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-2-0_gn_185001-201412.nc'
# --resolution  './TESTDIR_N96E_SSP585/qrparm.mask_n96e_eorca1_v2.2x_ESMF'
# -o './TESTDIR_JAN2019/Ndep_n96e_1893_1895' -n 'Ndep' --begin 1893 --end 1895
#
#
#
# These two commands produced ancillaries which were correct and had the same data (7/2061 looked at) as the original SSP126 ancillary created in December:
# python3 cmip6_ndep_jasmin_vn26_past_and_future.py -s './TESTDIR_N96E_SSP126/drynhx_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp126-1-0_gn_201501-209912.nc'
#   './TESTDIR_N96E_SSP126/drynoy_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp126-1-0_gn_201501-209912.nc'
#   './TESTDIR_N96E_SSP126/wetnhx_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp126-1-0_gn_201501-209912.nc'
#   './TESTDIR_N96E_SSP126/wetnoy_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp126-2-0_gn_201501-209912.nc'
#   --resolution  './TESTDIR_N96E_SSP126/qrparm.mask_n96e_eorca1_v2.2x_ESMF' -o './TESTDIR_JAN2019/Ndep_n96e_ssp126_2060_2062_appended'
#   -n 'Ndep' --begin 2060 --end 2062 -a true
#
# python3 cmip6_ndep_jasmin_vn26_past_and_future.py -s './TESTDIR_N96E_SSP126/drynhx_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp126-1-0_gn_201501-209912.nc'
#   './TESTDIR_N96E_SSP126/drynoy_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp126-1-0_gn_201501-209912.nc'
#   './TESTDIR_N96E_SSP126/wetnhx_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp126-1-0_gn_201501-209912.nc'
#   './TESTDIR_N96E_SSP126/wetnoy_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp126-2-0_gn_201501-209912.nc'
#   --resolution  './TESTDIR_N96E_SSP126/qrparm.mask_n96e_eorca1_v2.2x_ESMF' -o './TESTDIR_JAN2019/Ndep_n96e_ssp126_2060_2062'
#   -n 'Ndep' --begin 2060 --end 2062

# ---------------------------------------------------------------------------------
# Print intentions before processing
# ---------------------------------------------------------------------------------
def _process(start, end, name, sourcedir, output_file, resolution, appendheadtail):
    """State intentions before processing."""
    # Distinguish between timeslice (climatology for one year) / timeseries
    # based on whether an end argument is provided:
    if end is None:
        print(
            (
                "Process a timeslice of {} from files in directory {} for year {}, "
                "at {} resolution, writing result to {}."
                "NOTE: Likely only ever to be used for 1850 Preindustrial climatology:"
                "      If another year is requested the data from the CMIP"
                "      1850 climatology will be in the ancillary".format(
                    name, sourcedir, start, resolution, output_file
                )
            )
        )
    else:
        print(
            (
                "Process a timeseries of {} from files in directory {} between {} "
                "and {}, at resolution {} writing result to {}. AppendHeadTail = {}".format(
                    name, sourcedir, start, end, resolution, output_file, appendheadtail
                )
            )
        )


# -------------------------------------------------------------------------------------
# Procedure to selectively read in the climatology or timeseries files as appropriate
# -------------------------------------------------------------------------------------
def load_cube_clim_noclim(filenames, clim):
    """
    Discriminate between climatology and non-climatology data by whether the
    time coordinate has 12 points or not.
    """
    # Because data is loaded lazily, this is a sufficiently efficient
    # solution.
    cubes = iris.load(filenames)
    if clim:
        cubes = [cube for cube in cubes if len(cube.coord(axis="t").points) == 12]
    else:
        cubes = [cube for cube in cubes if len(cube.coord(axis="t").points) != 12]
    return cubes


# ---------------------------------------------------------------------------------------
# Main program
# ---------------------------------------------------------------------------------------
# def main():
if __name__ == "__main__":
    # __doc__ is the module docstring.
    arg_parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Add resolution as an argument
    arg_parser.add_argument(
        "-r",
        "--resolution",
        type=str,
        required=True,
        help="File containing grid to which the result "
        "will be gridded, and hence the definition "
        "of the output resolution.",
    )

    # Append another copy of the first year at the start and the last year at the end,
    # if necessary, to give the UM something to interpolate from/to when the
    # start/end date of the simulation matches that of the forcing dataset provided.
    arg_parser.add_argument(
        "-a",
        "--appendheadtail",
        type=str,
        default=None,
        help="Add second copy of start / end month at start / end of the ancillary",
    )

    arg_parser.add_argument(
        "--use-new-saver",
        action="store_true",
        required=False,
        help="If True, use the new saver.",
    )

    args = arg_parser.parse_args()
    _process(
        args.begin,
        args.end,
        args.name,
        args.sources,
        args.output,
        args.resolution,
        args.appendheadtail,
    )
    resolution = "n216e" if "n216e" in args.resolution else "n96e"
    appendHeadTail = args.appendheadtail
    gridFile = args.resolution
    if args.begin < 2014:
        pastOrFuture = "past"
        ancRefYear = "1850"
    else:
        pastOrFuture = "future"
        ancRefYear = str(args.begin)

    # If an end date has been given it's a timeseries, not a timeslice,
    # therefore not a climatology
    if args.end:
        clim = False
        type = "timeseries"

    # If no end date has been given, the required ancillary is a timeslice,
    # not a timeseries
    else:
        clim = True
        type = "timeslice"

    # Ancillary's attributes
    stash = "m01s00i447"
    description = "Nitrogen deposition field for driving UKESM1"


def main():

    # ---------------------------------
    # Create a TIMESERIES ancillary
    # ---------------------------------

    if type == "timeseries":
        tseriesCubeList = load_cube_clim_noclim(args.sources, 0)
        tseriesCube = (
            tseriesCubeList[0]
            + tseriesCubeList[1]
            + tseriesCubeList[2]
            + tseriesCubeList[3]
        )
        tseriesCube.coord("time").units = cf_units.Unit(
            "days since 1850-01-01 00:00:00", calendar="365_day"
        )
        if (pastOrFuture) == "past":
            tseriesCube.coord("time").units = cf_units.Unit(
                "days since 1850-01-01 00:00:00", calendar="360_day"
            )
            tseriesCube.coord("time").points = (
                tseriesCube.coord("time").points * 30.0
            )  # Needed for HIST / PI
        tseriesCube.coord("time").bounds = None

        with iris.FUTURE.context(cell_datetime_objects=True):
            constraint_ts = iris.Constraint(
                time=lambda cell: args.begin <= cell.point.year <= args.end
            )
            tseriesCube = tseriesCube.extract(constraint_ts)

        if (
            appendHeadTail
        ):  # Append start year to start year minus 1 and end year to end year plus 1, e.g.
            # for historical ancillary, copy 1850 field to 1849 and 2014 field to 2015, to
            # allow interpolatiion necessary for the UM to run from 1/1/1850 to 30/12/2014
            # This should only be TRUE if the entire timeseries of data provided is to be used,
            # e.g. the UM needs to be run from 1/1/1850 to 1/1/2015 and the data provided
            # covers the same period. Otherwise, simply create ancillaries which extend the
            # beyond the run period, e.g. for a run from 1950 to 1990, make ancillaries
            # from 1949 to 1991.

            with iris.FUTURE.context(cell_datetime_objects=True):
                constraint_start = iris.Constraint(
                    time=lambda cell: cell.point.year == args.begin + 1
                )
                constraint_end = iris.Constraint(
                    time=lambda cell: cell.point.year == args.end - 1
                )

                tseriesCube_start = tseriesCube.extract(constraint_start)
                tseriesCube_end = tseriesCube.extract(constraint_end)

            tseriesCube_startm1 = tseriesCube_start
            tseriesCube_endp1 = tseriesCube_end
            tseriesCube_startm1.coord("time").points = (
                tseriesCube_startm1.coord("time").points - 360.0
            )
            if tseriesCube_endp1 is not None:
                tseriesCube_endp1.coord("time").points = (
                    tseriesCube_endp1.coord("time").points + 360.0
                )
                tseriesCubes_startm1_to_endp1 = iris.cube.CubeList(
                    [
                        tseriesCube_startm1,
                        tseriesCube,
                        tseriesCube_endp1,
                    ]
                )
            else:
                # Need to pad two years at end
                with iris.FUTURE.context(cell_datetime_objects=True):
                    constraint_end = iris.Constraint(
                        time=lambda cell: cell.point.year == args.end - 2
                    )
                    tseriesCube_end = tseriesCube.extract(constraint_end)
                tseriesCube_endp1 = tseriesCube_end.copy()
                tseriesCube_endp2 = tseriesCube_end.copy()
                tseriesCube_endp1.coord("time").points = (
                    tseriesCube_endp1.coord("time").points + 720.0
                )
                tseriesCube_endp2.coord("time").points = (
                    tseriesCube_endp1.coord("time").points + 360.0
                )
                tseriesCubes_startm1_to_endp1 = iris.cube.CubeList(
                    [
                        tseriesCube_startm1,
                        tseriesCube,
                        tseriesCube_endp2,
                        tseriesCube_endp1,
                    ]
                )

            # Concatenate the three or four cubes: startm1, start --> end, endp1 or
            # startm1, start --> end, endp2, endp1
            tseriesCubes_startm1_to_endp1_cc = (
                tseriesCubes_startm1_to_endp1.concatenate_cube()[:]
            )
            tseriesCubes_startm1_to_endp1_cc.var_name = "n_dep"
            tseriesCubes_startm1_to_endp1_cc.long_name = (
                "NITROGEN DEPOSITION (kgN/m2/s)"
            )

            # New:
            tseriesCubes_startm1_to_endp1_cc.coord("time").points = [
                i * 30 + 15
                for i in range(
                    len(tseriesCubes_startm1_to_endp1_cc.coord("time").points)
                )
            ]

            # Regrid the source data to the UKESM grid
            beginm1 = str(args.begin - 1)
            tseriesCubeRG_out = regrid(tseriesCubes_startm1_to_endp1_cc)
            #           tseriesCubeRG_out.coord('time').units  = cf_units.Unit("days since 1850-01-01 00:00:00", calendar='360_day')
            tseriesCubeRG_out.coord("time").units = cf_units.Unit(
                "days since " + str(args.begin) + "-01-01 00:00:00", calendar="360_day"
            )

        else:

            tseriesCubeRG_out = regrid(tseriesCube)
            begin = str(args.begin)
            if (pastOrFuture) == "past":  # This worked for the 1871 1905 example
                tseriesCubeRG_out.coord("time").units = cf_units.Unit(
                    "days since " + ancRefYear + "-01-01 00:00:00", calendar="360_day"
                )
            else:  # Copied from append bit as was used for SSPs in December, so use here for SSPs
                # Easiest to use the 1t January of the starting year as the reference time, then build the points from  zero upwards.
                print(
                    "tseriesCubeRG_out.coord(time).points[0]  = ",
                    tseriesCubeRG_out.coord("time").points[0],
                )
                tseriesCubeRG_out.coord("time").units = cf_units.Unit(
                    "days since " + ancRefYear + "-01-01 00:00:00", calendar="360_day"
                )
                tseriesCubeRG_out.coord("time").points = [
                    i * 30 + 15
                    for i in range(len(tseriesCubeRG_out.coord("time").points))
                ]

        tseriesCubeRG_out.var_name = "n_dep"
        tseriesCubeRG_out.long_name = "NITROGEN DEPOSITION (kgN/m2/s)"
        tseriesCubeRG_out.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi(
            stash
        )
        tseriesCubeRG_out.attributes["grid_staggering"] = 6
        tseriesCubeRG_out.coord("time").bounds = None
        ants.utils.coord.guess_bounds(tseriesCubeRG_out.coord("time"), strict=True)
        tseriesCubeRG_out.data = np.clip(tseriesCubeRG_out.data, 0.0, float("Inf"))
        save_cube(tseriesCubeRG_out, args.output, args.use_new_saver)

        tseriesCubeRG_out.attributes["description"] = description
        print("Made " + args.output)

    # --------------------------------------------------------------------
    # Create timeslice ancillary from climatology
    # Currently only possible for preindustrial ancillary from 1850 from
    # climatology of monthly data
    # --------------------------------------------------------------------

    if type == "timeslice":

        tsliceCubeList = load_cube_clim_noclim(args.sources, clim)
        tsliceCube = (
            tsliceCubeList[0]
            + tsliceCubeList[1]
            + tsliceCubeList[2]
            + tsliceCubeList[3]
        )
        tsliceCube.var_name = "n_dep"
        tsliceCube.long_name = "NITROGEN DEPOSITION (kgN/m2/s)"
        tsliceCube.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi(stash)
        tsliceCube.attributes["grid_staggering"] = 6
        tsliceCube.coord("time").units = cf_units.Unit(
            "days since 1850-01-01 00:00:00", calendar="360_day"
        )
        tsliceCube.coord("time").points = tsliceCube.coord("time").points * 30.0
        tsliceCube.coord("time").bounds = None
        ants.utils.coord.guess_bounds(tsliceCube.coord("time"), strict=True)

        tsliceCubeRG = regrid(tsliceCube)
        tsliceCubeRG.data = np.clip(tsliceCubeRG.data, 0.0, float("Inf"))
        save_cube(tsliceCubeRG, args.output, args.use_new_saver)
        tsliceCubeRG.attributes["description"] = description
        print("Made " + args.output)


# -----------------------------------------------------------------------------------------------------
# Subroutine to regrid source data from its native grid to the model grid at the required resolution
# -----------------------------------------------------------------------------------------------------
def regrid(inCube):

    gridCube = ants.load_cube(gridFile, iris.AttributeConstraint(STASH="m01s00i030"))
    inCube.coord("latitude").coord_system = gridCube.coord("latitude").coord_system
    inCube.coord("longitude").coord_system = gridCube.coord("longitude").coord_system
    inCube.coord("longitude").circular = True
    gridCube.coord("longitude").circular = True
    ants.utils.cube.guess_horizontal_bounds(inCube)
    ants.utils.cube.guess_horizontal_bounds(gridCube)

    outCube = inCube.regrid(gridCube, iris.analysis.AreaWeighted())
    return outCube


def save_cube(cube, filename, use_new_saver):
    if use_new_saver:
        save.ancil(cube, filename)
        save.netcdf(cube, filename)
    else:
        ants.save(cube, filename, saver="ancil")
    print("Made " + filename)


if __name__ == "__main__":
    main()
