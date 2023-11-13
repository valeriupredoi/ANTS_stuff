"""
Create the ancillaries needed to represent nitrogen deposition in UKESM1 for the CMIP6/DEC historical simulation 
00447    0 NITROGEN DEPOSITION (kgN/m2/s)
"""
import datetime  # SKL added

import cf_units
import iris
import iris.analysis
import iris.fileformats
import numpy as np
from iris.time import PartialDateTime  # SKL added

import ants
import ants.decomposition as decomp
import ants.fileformats

iris.FUTURE.netcdf_promote = True  # SKL added, got rid of a warning message

# Read in the resolution-specific land-sea mask

n96eFile = "/home/users/sliddicoat/CMIP6_FORCING/RAW_DATA/GRIDFILES/atmos/n96e_eorca1_v2.2x_ESMF/land_sea_mask/etop01/qrparm.mask_n96e_eorca1_v2.2x_ESMF"
n216eFile = "/group_workspaces/jasmin2/tids/UM/ancil/atmos/n216e/orca025/land_sea_mask/etop01/v2/qrparm.mask"

gridFiles = {}
gridFiles["n96e"] = n96eFile
gridFiles["n216e"] = n216eFile

resolution = "n96e"  # resolution of ancillary to create
gridFile = gridFiles[resolution]

print("gridFile = ", gridFile)


# Create output filenames

stash = "m01s00i447"

# Timeslices

pi_ancil_name = (
    "/home/users/sliddicoat/CMIP6_FORCING/ANCILS/Ndep/ndep_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-2-0_gn_185001-185012-clim."
    + resolution
    + ".anc"
)
hist_ancil_name = (
    "/home/users/sliddicoat/CMIP6_FORCING/ANCILS/Ndep/ndep_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-2-0_gn_184901-201512."
    + resolution
    + ".anc"
)
ts2014_ancil_name = (
    "/home/users/sliddicoat/CMIP6_FORCING/ANCILS/Ndep/ndep_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-2-0_gn_201401-201412."
    + resolution
    + ".anc"
)

piSwitch = True  # Preindustrial climatology - uses separate files from the Historical
histSwitch = True  # Historical timeseries
histAppendHeadTail = True  # If true, copy 1850 data to 1849 and 2014 data to 2015 to allow UM to have data to interpolate to at start / end of run

tsSwitch = True  # timeslice using data from the historical files
tsYearStart = 2014  # timeslice - start year
tsYearEnd = 2014  # timeslice - end year


def main():

    # Quick test to make sure all cubes are loading having changed from one directory in my space to the 4 individual dry/wet/nhx/noy directories in the inputs4MIPs mirror on Jasmin

    Files = "*_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-2-0_gn_185001-185012-clim.nc"

    inDir1 = "/home/users/sliddicoat/CMIP6_FORCING/RAW_DATA/Ndep/"
    testList1 = ants.load(inDir1 + Files)
    print()
    print("testList1 = ", testList1)

    inDir2 = "/group_workspaces/jasmin2/tids/CMIP6_ANCIL/data/inputs4MIPs_2017-06-08/UReading/surfaceFluxes/CMIP/NCAR-CCMI-2-0/mon/*/gn/v20161207/"
    testList2 = ants.load(inDir2 + Files)
    print()
    print("testList2 = ", testList2)

    # testList1 =  0: drynhx / (kg m-2 s-1)               (time: 12; latitude: 96; longitude: 144)
    # 1: wetnhx / (kg m-2 s-1)               (time: 12; latitude: 96; longitude: 144)
    # 2: tendency_of_atmosphere_mass_content_of_noy_expressed_as_nitrogen_due_to_dry_deposition / (kg m-2 s-1) (time: 12; latitude: 96; longitude: 144)
    # 3: tendency_of_atmosphere_mass_content_of_noy_expressed_as_nitrogen_due_to_wet_deposition / (kg m-2 s-1) (time: 12; latitude: 96; longitude: 144)

    # testList2 =  0: drynhx / (kg m-2 s-1)               (time: 12; latitude: 96; longitude: 144)
    # 1: wetnhx / (kg m-2 s-1)               (time: 12; latitude: 96; longitude: 144)
    # 2: tendency_of_atmosphere_mass_content_of_noy_expressed_as_nitrogen_due_to_dry_deposition / (kg m-2 s-1) (time: 12; latitude: 96; longitude: 144)
    # 3: tendency_of_atmosphere_mass_content_of_noy_expressed_as_nitrogen_due_to_wet_deposition / (kg m-2 s-1) (time: 12; latitude: 96; longitude: 144)

    inDirs = "/group_workspaces/jasmin2/tids/CMIP6_ANCIL/data/inputs4MIPs_2017-06-08/UReading/surfaceFluxes/CMIP/NCAR-CCMI-2-0/mon/*/gn/v20161207/"
    piFiles = "*_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-2-0_gn_185001-185012-clim.nc"
    histFiles = "*_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-2-0_gn_185001-201412.nc"

    if piSwitch:
        piCubeList = ants.load(inDirs + piFiles)
        print()
        print("piCubeList ", piCubeList)
        print()
        piCube = piCubeList[0] + piCubeList[1] + piCubeList[2] + piCubeList[3]
        piCube = regrid(piCube)
        print(piCube)
        piCube.var_name = "n_dep"
        piCube.long_name = "NITROGEN DEPOSITION (kgN/m2/s)"
        piCube.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi(stash)
        piCube.attributes["grid_staggering"] = 6
        print()
        print("piCube.attributes = ", piCube.attributes)
        print()
        piCube.coord("time").units = cf_units.Unit(
            "days since 1850-01-01 00:00:00", calendar="360_day"
        )
        print()
        print("piCube.coord(time).units =", piCube.coord("time").units)
        print()
        piCube.coord("time").points = piCube.coord("time").points * 30.0
        print()
        print("piCube.coord(time).points =", piCube.coord("time").points)
        print()
        piCube.coord("time").bounds = None
        ants.utils.coord.guess_bounds(piCube.coord("time"), strict=True)  # SKL added
        print(piCube)
        print(piCube.__class__.__name__)
        piCube.data = np.clip(piCube.data, 0.0, float("Inf"))
        ants.save(piCube, pi_ancil_name, saver="ancil")
        print("made " + pi_ancil_name)

    if histSwitch:
        histCubeList = ants.load(inDirs + histFiles)
        print()
        print("histCubeList = ", histCubeList)
        print()
        histCube = histCubeList[0] + histCubeList[1] + histCubeList[2] + histCubeList[3]

        histCube.coord("time").units = cf_units.Unit(
            "days since 1850-01-01 00:00:00", calendar="360_day"
        )

        histCube.coord("time").points = (
            histCube.coord("time").points * 30.0
        )  # Necessary as otherwise dates
        # go up daily from 1850/01/01 to 1855 rather than monthly to 2014
        histCube.coord("time").bounds = None

        if histAppendHeadTail:

            with iris.FUTURE.context(cell_datetime_objects=True):
                constraint_1850 = iris.Constraint(
                    time=lambda cell: 1850 <= cell.point.year <= 1850
                )
                constraint_2014 = iris.Constraint(
                    time=lambda cell: 2014 <= cell.point.year <= 2014
                )
                histCube_1850 = histCube.extract(constraint_1850)
                histCube_2014 = histCube.extract(constraint_2014)

            # Copy 2014 cube to 2015 and 1850 cube to 1849

            histCube_1849 = histCube_1850
            histCube_2015 = histCube_2014

            histCube_1849.coord("time").points = (
                histCube_1849.coord("time").points - 360.0
            )
            histCube_2015.coord("time").points = (
                histCube_2015.coord("time").points + 360.0
            )
            print()
            print("histCube_1850.coord(time) = ", histCube_1850.coord("time"))
            print("histCube_2014.coord(time) = ", histCube_2014.coord("time"))
            print()
            print()
            print(
                "histCube_1850.coord(time).points = ",
                histCube_1850.coord("time").points,
            )
            print(
                "histCube_2014.coord(time).points = ",
                histCube_2014.coord("time").points,
            )
            print()
            print()

            # Concatenate the three cubes: 1849, 1850-2014, 2015

            histCubes_1849_to_2015 = iris.cube.CubeList(
                [histCube_1849, histCube, histCube_2015]
            )
            histCubes_1849_to_2015_cc = histCubes_1849_to_2015.concatenate_cube()[:]

            print()
            print("histCubes_1849_to_2015                = ", histCubes_1849_to_2015)
            print("histCubes_1849_to_2015_cc             = ", histCubes_1849_to_2015_cc)
            print()
            print(
                "histCubes_1849_to_2015_cc.coord(time) = ",
                histCubes_1849_to_2015_cc.coord("time"),
            )
            print(
                "histCubes_1849_to_2015_cc.coord(time).points = ",
                histCubes_1849_to_2015_cc.coord("time").points,
            )
            print()

            # Regrid to the target grid

            histCube_out = regrid(histCubes_1849_to_2015_cc)

        else:
            histCube_out = regrid(histCube)

    histCube_out.var_name = "n_dep"
    histCube_out.long_name = "NITROGEN DEPOSITION (kgN/m2/s)"
    histCube_out.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi(stash)
    histCube_out.attributes["grid_staggering"] = 6

    print()
    print("histCube_out.coord(time)  = ", histCube_out.coord("time"))
    print()

    histCube_out.coord("time").bounds = None
    ants.utils.coord.guess_bounds(histCube_out.coord("time"), strict=True)

    print()
    print(
        "histCube_out.coord(time).points after guess bounds = ",
        histCube_out.coord("time").points,
    )
    print()
    print(histCube)
    print(histCube.__class__.__name__)
    print()

    histCube.data = np.clip(histCube.data, 0.0, float("Inf"))
    ants.save(histCube_out, hist_ancil_name, saver="ancil")
    print("made " + hist_ancil_name)


def regrid(inCube):  # ,gridCube):
    #    """
    #    regrid the data onto the HG3 n96e grid
    #    new grid is specified by "gridFile"
    #    gridFile is a fractional area ancillary
    #    gridFile is not passed into the routine because the domain decomposition can't handle cubes with different shapes
    #    """
    gridCube = ants.load_cube(gridFile, iris.AttributeConstraint(STASH="m01s00i030"))
    inCube.coord("latitude").coord_system = gridCube.coord("latitude").coord_system
    inCube.coord("longitude").coord_system = gridCube.coord("longitude").coord_system
    inCube.coord("longitude").circular = True
    gridCube.coord("longitude").circular = True
    ants.utils.cube.guess_horizontal_bounds(inCube)  # SKL added
    ants.utils.cube.guess_horizontal_bounds(gridCube)  # SKL added

    outCube = inCube.regrid(gridCube, iris.analysis.AreaWeighted())
    return outCube


if __name__ == "__main__":
    main()
