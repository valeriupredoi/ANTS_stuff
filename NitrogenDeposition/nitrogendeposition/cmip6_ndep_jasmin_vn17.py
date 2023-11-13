"""
Create the ancillaries needed to represent nitrogen deposition in UKESM1 for the CMIP6/DEC historical simulation 
00447    0 NITROGEN DEPOSITION (kgN/m2/s)
"""
import datetime

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

iris.FUTURE.netcdf_promote = True

a4port = (8.27, 11.69)  # A4 paper portrait in inches
a4land = (11.69, 8.27)  # A4 paper landscape in inches
colour = ["green", "blue", "orange", "red", "purple", "black"]

# Read in the resolution-specific land-sea mask

n96eFile = "/home/users/sliddicoat/CMIP6_FORCING/RAW_DATA/GRIDFILES/atmos/n96e_eorca1_v2.2x_ESMF/land_sea_mask/etop01/qrparm.mask_n96e_eorca1_v2.2x_ESMF"
n216eFile = "/group_workspaces/jasmin2/tids/UM/ancil/atmos/n216e/orca025/land_sea_mask/etop01/v2/qrparm.mask"

gridFiles = {}
gridFiles["n96e"] = n96eFile
gridFiles["n216e"] = n216eFile

resolution = "n96e"  # resolution of ancillary to create
resolution = "n216e"  # resolution of ancillary to create

gridFile = gridFiles[resolution]

# Specify ancillaries to produce and set output filenames and attributes

ancdir = "/group_workspaces/jasmin2/tids/CMIP6_ANCIL/users/sliddicoat/ANCILS/Ndep/"

piSwitch = True  # Preindustrial climatology - uses separate files from the Historical
pi_ancil_name = (
    ancdir
    + "ndep_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-2-0_gn_185001-185012-clim."
    + resolution
    + ".anc"
)

histSwitch = True  # Historical timeseries
histAppendHeadTail = True  # If true, copy 1850 data to 1849 and 2014 data to 2015 to allow UM to have data to interpolate to at start / end of run
hist_ancil_name = (
    ancdir
    + "ndep_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-2-0_gn_184901-201512."
    + resolution
    + ".anc"
)

tsSwitch = True  # timeslice using data from the historical files
tsStartYear = 2014  # timeslice - start year
tsEndYear = 2014  # timeslice - end year
if tsStartYear == tsEndYear:
    ts_ancil_name = (
        ancdir
        + "ndep_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-2-0_gn_ts_"
        + str(tsStartYear)
        + "."
        + resolution
        + ".anc"
    )
else:
    ts_ancil_name = (
        ancdir
        + "ndep_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-2-0_gn_ts_"
        + str(tsStartYear)
        + "_"
        + str(tsEndYear)
        + "."
        + resolution
        + ".anc"
    )

stash = "m01s00i447"
description = "Nitrogen deposition field for driving UKESM1"


def main():

    inDirs = "/group_workspaces/jasmin2/tids/CMIP6_ANCIL/data/inputs4MIPs_2017-06-08/UReading/surfaceFluxes/CMIP/NCAR-CCMI-2-0/mon/*/gn/v20161207/"
    piFiles = "*_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-2-0_gn_185001-185012-clim.nc"
    histFiles = "*_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-2-0_gn_185001-201412.nc"

    # ----------------------------------------------------------------
    # Create preindustrial ancillary from the 1850 climatology files
    # ----------------------------------------------------------------

    if piSwitch:
        piCubeList = ants.load(inDirs + piFiles)
        piCube = piCubeList[0] + piCubeList[1] + piCubeList[2] + piCubeList[3]
        piCubeRG = regrid(piCube)
        piCubeRG.var_name = "n_dep"
        piCubeRG.long_name = "NITROGEN DEPOSITION (kgN/m2/s)"
        piCubeRG.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi(stash)
        piCubeRG.attributes["grid_staggering"] = 6
        piCubeRG.coord("time").units = cf_units.Unit(
            "days since 1850-01-01 00:00:00", calendar="360_day"
        )
        piCubeRG.coord("time").points = piCubeRG.coord("time").points * 30.0
        piCubeRG.coord("time").bounds = None
        ants.utils.coord.guess_bounds(piCubeRG.coord("time"), strict=True)
        piCubeRG.data = np.clip(piCubeRG.data, 0.0, float("Inf"))
        ants.save(piCubeRG, pi_ancil_name, saver="ancil")
        print("made " + pi_ancil_name)
        piCubeRG.attributes["description"] = description
        ants.fileformats.netcdf.save(piCubeRG, pi_ancil_name + ".nc")
        print("made " + pi_ancil_name + ".nc")

    # ---------------------------------
    # Create the historical ancillary
    # ---------------------------------

    if histSwitch:
        histCubeList = ants.load(inDirs + histFiles)
        histCube = histCubeList[0] + histCubeList[1] + histCubeList[2] + histCubeList[3]

        histCube.coord("time").units = cf_units.Unit(
            "days since 1850-01-01 00:00:00", calendar="360_day"
        )
        histCube.coord("time").points = histCube.coord("time").points * 30.0
        histCube.coord("time").bounds = None

        if (
            histAppendHeadTail
        ):  # Copy 1850 field to 1849 and 2014 field to 2015 allow interpolatiion necessary
            # for the UM to run from 1/1/1850 to 30/12/2014

            with iris.FUTURE.context(cell_datetime_objects=True):
                constraint_1850 = iris.Constraint(
                    time=lambda cell: 1850 <= cell.point.year <= 1850
                )
                constraint_2014 = iris.Constraint(
                    time=lambda cell: 2014 <= cell.point.year <= 2014
                )
                histCube_1850 = histCube.extract(constraint_1850)
                histCube_2014 = histCube.extract(constraint_2014)

            histCube_1849 = histCube_1850
            histCube_2015 = histCube_2014

            histCube_1849.coord("time").points = (
                histCube_1849.coord("time").points - 360.0
            )
            histCube_2015.coord("time").points = (
                histCube_2015.coord("time").points + 360.0
            )

            # Concatenate the three cubes: 1849, 1850-2014, 2015

            histCubes_1849_to_2015 = iris.cube.CubeList(
                [histCube_1849, histCube, histCube_2015]
            )
            histCubes_1849_to_2015_cc = histCubes_1849_to_2015.concatenate_cube()[:]
            histCubes_1849_to_2015_cc.var_name = "n_dep"

            # Regrid the source data to the UKESM grid

            histCubeRG_out = regrid(histCubes_1849_to_2015_cc)

        else:
            histCubeRG_out = regrid(histCube)

    histCubeRG_out.coord("time").units = cf_units.Unit(
        "days since 1850-01-01 00:00:00", calendar="360_day"
    )
    histCubeRG_out.var_name = "n_dep"
    histCubeRG_out.long_name = "NITROGEN DEPOSITION (kgN/m2/s)"
    histCubeRG_out.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi(stash)
    histCubeRG_out.attributes["grid_staggering"] = 6
    histCubeRG_out.coord("time").bounds = None
    ants.utils.coord.guess_bounds(histCubeRG_out.coord("time"), strict=True)
    histCubeRG_out.data = np.clip(histCubeRG_out.data, 0.0, float("Inf"))
    ants.save(histCubeRG_out, hist_ancil_name, saver="ancil")
    print("made " + hist_ancil_name)

    histCubeRG_out.attributes["description"] = description
    ants.fileformats.netcdf.save(histCubeRG_out, hist_ancil_name + ".nc")
    print("made " + hist_ancil_name + ".nc")

    histCubes_1849_to_2015_cc.long_name = "NITROGEN DEPOSITION (kgN/m2/s)"

    #  Produce maps and timeseries of the original and regridded data
    plot_ndep_anc(
        histCubes_1849_to_2015_cc[:], histCubeRG_out, ["Input4MIPs", "Ancil."]
    )

    # -----------------------------
    # Create timeslice ancillary
    # -----------------------------

    if tsSwitch:
        tsCubeList = ants.load(inDirs + histFiles)
        tsCube = tsCubeList[0] + tsCubeList[1] + tsCubeList[2] + tsCubeList[3]
        tsCube.var_name = "n_dep"
        tsCube.long_name = "NITROGEN DEPOSITION (kgN/m2/s)"
        tsCube.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi(stash)
        tsCube.attributes["grid_staggering"] = 6
        tsCube.coord("time").units = cf_units.Unit(
            "days since 1850-01-01 00:00:00", calendar="360_day"
        )
        tsCube.coord("time").points = tsCube.coord("time").points * 30.0
        tsCube.coord("time").bounds = None
        ants.utils.coord.guess_bounds(tsCube.coord("time"), strict=True)

        with iris.FUTURE.context(cell_datetime_objects=True):
            constraint_ts = iris.Constraint(
                time=lambda cell: tsStartYear <= cell.point.year <= tsEndYear
            )
            tsCube = tsCube.extract(constraint_ts)

        tsCubeRG = regrid(tsCube)
        tsCubeRG.data = np.clip(tsCubeRG.data, 0.0, float("Inf"))
        ants.save(tsCubeRG, ts_ancil_name, saver="ancil")
        print("made " + ts_ancil_name)
        tsCubeRG.attributes["description"] = description
        ants.fileformats.netcdf.save(tsCubeRG, ts_ancil_name + ".nc")
        print("made " + ts_ancil_name + ".nc")


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

    fig.savefig(ancdir + "maps_Ndep_" + resolution + ".png")

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

    fig.savefig(ancdir + "timeseries_Ndep_" + resolution + ".png")
    iplt.show()


if __name__ == "__main__":
    main()
