#!/usr/bin/env python
"""
Create the ancillaries needed to represent nitrogen deposition in UKESM1 for the CMIP6/DEC historical simulation 
00447    0 NITROGEN DEPOSITION (kgN/m2/s)
"""
import datetime
import glob

import ants
import ants.decomposition as decomp
import ants.fileformats
import ants.io.save as save
import cf_units

import argparse

import iris
import iris.analysis
import iris.fileformats
import iris.plot as iplt
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np
import pylab
from iris.time import PartialDateTime

print ("Here")
iris.FUTURE.netcdf_promote = True

a4port = (8.27, 11.69)  # A4 paper portrait in inches
a4land = (11.69, 8.27)  # A4 paper landscape in inches
colour = ["green", "blue", "orange", "red", "purple", "black"]

# SKL 7th November 2018
#
# New version of the script (vn18) to make use of the command line interface script, cli.py (imported above),
# to allow this script to be run from the command line. Example calls below and the output showing the
# resulting files created.
#
# $ python2.7 cmip6_ndep_jasmin_vn21.py -d './TESTDIR_N96E/' -r 'n96e' -o './TESTDIR_N96E/n96e_ancillary_1914_1918' -n 'Ndep' --begin 1914 --end 1918
# Made ./TESTDIR_N96E/n96e_ancillary_1914_1918.anc
# Made ./TESTDIR_N96E/n96e_ancillary_1914_1918.nc
#
# $ python2.7 cmip6_ndep_jasmin_vn21.py -d './TESTDIR_N96E/' -r 'n96e' -o './TESTDIR_N96E/n96e_ancillary_1850_clim' -n 'Ndep' --begin 1850
# Made ./TESTDIR_N96E/n96e_ancillary_1850_clim.anc
# Made ./TESTDIR_N96E/n96e_ancillary_1850_clim.nc
#
# $ python2.7 cmip6_ndep_jasmin_vn21.py -d './TESTDIR_N96E/' -r 'n96e' -o './TESTDIR_N96E/n96e_ancillary_1849_2015' -n 'Ndep' --begin 1850 --end 2014 --appendheadtail 1
# Made ./TESTDIR_N96E/n96e_ancillary_1849_2015.anc
# Made ./TESTDIR_N96E/n96e_ancillary_1849_2015.nc
#
# $ python2.7 cmip6_ndep_jasmin_vn21.py -d './TESTDIR_N96E/' -r 'n96e' -o './TESTDIR_N96E/n96e_ancillary_1850_2014' -n 'Ndep' --begin 1850 --end 2014
# Made ./TESTDIR_N96E/n96e_ancillary_1850_2014.anc
# Made ./TESTDIR_N96E/n96e_ancillary_1850_2014.nc
#
# python2.7 cmip6_ndep_jasmin_vn22_ssp585.py -s '/TESTDIR_N96E_SSP585/drynhx_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp585-1-0_gn_201501-209912.nc' './TESTDIR_N96E_SSP585/drynoy_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp585-1-0_gn_201501-209912.nc' './TESTDIR_N96E_SSP585/wetnhx_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp585-1-0_gn_201501-209912.nc' './TESTDIR_N96E_SSP585/wetnoy_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp585-2-0_gn_201501-209912.nc' --resolution  './TESTDIR_N96E_SSP585/qrparm.mask_n96e_eorca1_v2.2x_ESMF' -o './TESTDIR_N96E_SSP585/n96e_ancillary_2015_2100_v2' -n 'Ndep' --begin 2015 --end 2100 -a true
#

# ---------------------------------------------------------------------------------
# Print intentions before processing
# ---------------------------------------------------------------------------------
def _process(start, end, name, sourcedir, output_file, resolution, appendheadtail):
    """State intentions before processing."""
    # Distinguish between timeslice (climatology for one year) / timeseries
    # based on whether an end argument is provided:
    if end is None:
        print((
            "Process a timeslice of {} from files in directory {} for year {}, "
            "at {} resolution, writing result to {}."
            "NOTE: Likely only ever to be used for 1850 Preindustrial climatology:"
            "      If another year is requested the data from the CMIP"
            "      1850 climatology will be in the ancillary".format(
                name, sourcedir, start, resolution, output_file
            )
        ))
    else:
        print((
            "Process a timeseries of {} from files in directory {} between {} "
            "and {}, at resolution {} writing result to {}. AppendHeadTail = {}".format(
                name, sourcedir, start, end, resolution, output_file, appendheadtail
            )
        ))


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
    print("IN LOAD_CUBE_CLIM_NOCLIM: filenames = ", filenames)
    print("IN LOAD_CUBE_CLIM_NOCLIM: cubes     = ", cubes)
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
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

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

    # Append another copy of the first year at the start and the last year at the end
    # if necessary. Probably only reqiured for Historical as CMIP dataset
    # goes from 1850 to 2014: for an experiment from 1850 to 2014, the UM requires
    # data from the end of 1849 to the start of 2015 for time interpolation.
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

    print("clim = ", clim)
    print("gridFile = ", gridFile)

    # Specify ancillaries to produce and set output filenames and attributes

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
        print("tseriesCube AAA = ", tseriesCube)
        print("tseriesCube.coord(time).points AAA = ", tseriesCube.coord("time").points)
        print("len(tseriesCube.coord(time).points) AAA = ", len(
            tseriesCube.coord("time").points
        ))
        #        tseriesCube.coord('time').units  = cf_units.Unit("days since 1850-01-01 00:00:00", calendar='360_day')
        tseriesCube.coord("time").units = cf_units.Unit(
            "days since 1850-01-01 00:00:00", calendar="365_day"
        )
        #  tseriesCube.coord('time').points = tseriesCube.coord('time').points*30.0	# Needed for HIST / PI
        tseriesCube.coord("time").bounds = None

        with iris.FUTURE.context(cell_datetime_objects=True):
            #           constraint_ts = iris.Constraint(time=lambda cell: args.begin <= cell.point.year <= args.begin)
            constraint_ts = iris.Constraint(
                time=lambda cell: args.begin <= cell.point.year <= args.end
            )
            print("args.begin = ", args.begin)
            print("args.end   = ", args.end)
            print("args.sources = ", args.sources)
            tseriesCube = tseriesCube.extract(constraint_ts)
            print("tseriesCube BBB = ", tseriesCube)
            print("tseriesCube.coord(time).points BBB = ", tseriesCube.coord(
                "time"
            ).points)
            print("len(tseriesCube.coord(time).points) BBB = ", len(
                tseriesCube.coord("time").points
            ))
            #          tseriesCube.coord('time').points =[ i*30 +15 for i in range(len(tseriesCube.coord('time').points)) ]
            print("tseriesCube.coord(time).points CCC = ", tseriesCube.coord(
                "time"
            ).points)
            print("len(tseriesCube.coord(time).points) CCC = ", len(
                tseriesCube.coord("time").points
            ))

        if (
            appendHeadTail
        ):  # Append start year to start year minus 1 and end year to end year plus 1, e.g.
            # for historical ancillary, copy 1850 field to 1849 and 2014 field to 2015, to
            # allow interpolatiion necessary for the UM to run from 1/1/1850 to 30/12/2014

            with iris.FUTURE.context(cell_datetime_objects=True):
                constraint_start = iris.Constraint(
                    time=lambda cell: args.begin <= cell.point.year <= args.begin
                )
                constraint_end = iris.Constraint(
                    time=lambda cell: args.end <= cell.point.year <= args.end
                )
                print("args.begin, args.end  = ", args.begin, args.end)
                tseriesCube_start = tseriesCube.extract(constraint_start)
                tseriesCube_end = tseriesCube.extract(constraint_end)
                print("tseriesCube_start = ", tseriesCube_start)
                print("tseriesCube_end = ", tseriesCube_end)
            #     print 'tseriesCube_stat = ', tseriesCube_start
            tseriesCube_startm1 = tseriesCube_start
            tseriesCube_endp1 = tseriesCube_end
            #           print 'DOING STARTM1: tseriesCube_startm1.coord(time).points = ',tseriesCube_startm1.coord('time').points
            tseriesCube_startm1.coord("time").points = (
                tseriesCube_startm1.coord("time").points - 360.0
            )
            tseriesCube_endp1.coord("time").points = (
                tseriesCube_endp1.coord("time").points + 360.0
            )

            # Concatenate the three cubes: startm1, start --> end, endp1

            tseriesCubes_startm1_to_endp1 = iris.cube.CubeList(
                [tseriesCube_startm1, tseriesCube, tseriesCube_endp1]
            )
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
            print("beginm11 = ", beginm1)
            tseriesCubeRG_out = regrid(tseriesCubes_startm1_to_endp1_cc)
            tseriesCubeRG_out.coord("time").units = cf_units.Unit(
                "days since " + beginm1 + "-01-01 00:00:00", calendar="360_day"
            )

        else:
            print("tseriesCube DDD = ", tseriesCube)
            tseriesCubeRG_out = regrid(tseriesCube)

            #        tseriesCube.coord('time').points =[ i*30 +15 for i in range(len(tseriesCube.coord('time').points)) ]
            #        tseriesCube.coord('time').points =[ i*30 +15 for i in range(len(tseriesCube.coord('time').points)) ]
            #        tseriesCubeRG_out.coord('time').units  = cf_units.Unit("days since 1850-01-01 00:00:00", calendar='360_day')
            begin = str(args.begin)
            tseriesCubeRG_out.coord("time").units = cf_units.Unit(
                "days since " + begin + "-01-01 00:00:00", calendar="360_day"
            )

        tseriesCubeRG_out.var_name = "n_dep"
        tseriesCubeRG_out.long_name = "NITROGEN DEPOSITION (kgN/m2/s)"
        tseriesCubeRG_out.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi(
            stash
        )
        tseriesCubeRG_out.attributes["grid_staggering"] = 6
        tseriesCubeRG_out.coord("time").bounds = None
        ants.utils.coord.guess_bounds(tseriesCubeRG_out.coord("time"), strict=True)
        tseriesCubeRG_out.data = np.clip(tseriesCubeRG_out.data, 0.0, float("Inf"))
        save_cube(tseriesCubeRG_out, args.output + ".anc", args.use_new_saver)

        tseriesCubeRG_out.attributes["description"] = description
        print("Made " + args.output + ".nc")

        #  Produce maps and timeseries of the original and regridded data
        if appendHeadTail:
            plot_ndep_anc(
                tseriesCubes_startm1_to_endp1_cc[:],
                tseriesCubeRG_out,
                ["Input4MIPs", "Ancil."],
            )
        else:
            plot_ndep_anc(tseriesCube[:], tseriesCubeRG_out, ["Input4MIPs", "Ancil."])

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
        save_cube(tsliceCubeRG, args.output + ".anc", args.use_new_saver)
        tsliceCubeRG.attributes["description"] = description
        print("Made " + args.output + ".nc")
        # plot_ndep_anc(tsliceCube[:],tsliceCubeRG,["Input4MIPs","Ancil."])


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


# ------------------------------------------------------------------------------------
# Plot maps and timeseries of source data and regridded data for comparison
# ------------------------------------------------------------------------------------
def plot_ndep_anc(cl1, cl2, names):

    fig = plt.figure(figsize=a4land)

    # Plot maps from start of timeseries
    ax1 = fig.add_subplot(2, 2, 1)
    qplt.contourf(cl1[0, :, :])

    ax2 = fig.add_subplot(2, 2, 2)
    qplt.contourf(cl2[0, :, :])

    # Plot maps from end of timeseries
    ax1 = fig.add_subplot(2, 2, 3)
    qplt.contourf(cl1[-1, :, :])

    ax2 = fig.add_subplot(2, 2, 4)
    qplt.contourf(cl2[-1, :, :])

    fig.savefig(args.output + "_maps_Ndep.png")

    # Plot timeseries of global means
    fig = plt.figure(figsize=a4land)
    ants.utils.cube.guess_horizontal_bounds(cl1)
    ants.utils.cube.guess_horizontal_bounds(cl2)
    cl1glob = []
    cl2glob = []

    cl1glob.append(
        cl1.collapsed(
            ["longitude", "latitude"],
            iris.analysis.SUM,
            weights=iris.analysis.cartography.area_weights(cl1),
        )
    )

    cl1glob[0].data = cl1glob[0].data * 1.0e-12 * 60.0 * 60.0 * 24.0 * 360.0

    cl2glob.append(
        cl2.collapsed(
            ["longitude", "latitude"],
            iris.analysis.SUM,
            weights=iris.analysis.cartography.area_weights(cl2),
        )
    )

    cl2glob[0].data = cl2glob[0].data * 1.0e-12 * 60.0 * 60.0 * 24.0 * 360.0
    cl1glob[0].coord("time").units = cf_units.Unit(
        "days since 1850-01-01 00:00:00", calendar="360_day"
    )
    cl2glob[0].coord("time").units = cf_units.Unit(
        "days since 1850-01-01 00:00:00", calendar="360_day"
    )

    qp = qplt.plot(cl1glob[0], label="Input4MIPs raw data", color=colour[1])
    qp = qplt.plot(
        cl2glob[0],
        label="Regridded onto " + resolution + " grid",
        color=colour[2],
        linestyle="--",
    )

    pylab.xlim(1850, 2050)
    plt.axis("tight")

    pylab.ylabel("Total Nitrogen deposition [GtN yr-$1$]")
    plt.legend(loc="upper left", prop={"size": 10})

    fig.savefig(args.output + "_timeseries_Ndep.png")
    iplt.show()


def save_cube(cube, filename, use_new_saver):
    if use_new_saver:
        save.ancil(cube, filename)
        save.netcdf(cube, filename)
    else:
        ants.save(cube, filename, saver="ancil")
    print("Made " + filename)


if __name__ == "__main__":
    main()
