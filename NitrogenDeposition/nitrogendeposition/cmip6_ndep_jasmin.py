"""
Create the ancillaries needed to represent nitrogen deposition in UKESM1 for the CMIP6/DEC historical simulation 
00447    0 NITROGEN DEPOSITION (kgN/m2/s)
"""
import ants
import ants.decomposition as decomp
import ants.fileformats
import cf_units
import iris
import iris.analysis
import iris.fileformats
import numpy as np

iris.FUTURE.netcdf_promote = True  # SKL added, got rid of an warning message

n96eFile = "/home/users/sliddicoat/CMIP6_FORCING/RAW_DATA/GRIDFILES/um1/ancil/atmos/n96e/orca1/land_sea_mask/etop01/v0/qrparm.mask"
n216eFile = "/home/users/sliddicoat/CMIP6_FORCING/RAW_DATA/GRIDFILES/um1/ancil/atmos/n216e/orca1/land_sea_mask/etop01/v0/qrparm.mask"
n512eFile = "/home/users/sliddicoat/CMIP6_FORCING/RAW_DATA/GRIDFILES/um1/ancil/atmos/n5126e/orca1/land_sea_mask/etop01/v0/qrparm.mask"

gridFile = n96eFile

stash = "m01s00i447"
pi_ancil_name = "/home/users/sliddicoat/CMIP6_FORCING/ANCILS/Ndep/ndep_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-1-0_gr_185001-185012.n96.anc"
hist_ancil_name = "/home/users/sliddicoat/CMIP6_FORCING/ANCILS/Ndep/ndep_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-1-0_gr_185001-201412.n96.anc"
piSwitch = True
histSwitch = True


# help(ants)


def main():
    #     laiCube = ants.load_cube("/hpc/projects/um1/ancil/atmos/n96e/orca1/vegetation/func_type_modis/v2/qrparm.veg.func",
    #                              iris.AttributeConstraint(STASH='m01s00i217'))
    #     print laiCube
    #     print laiCube.coord("time").units
    #     print laiCube.coord("time").points
    #     print laiCube.coord("time").coord_system
    #     exit()

    #    inDir = "/data/users/hadsl/UKESM1/NITROGEN_DEPOSITION_DATASETS/"
    inDir = "/home/users/sliddicoat/CMIP6_FORCING/RAW_DATA/Ndep/"
    piFiles = "*_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-1-0_gr_185001-185012.nc"
    histFiles = "*_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-1-0_gr_185001-201412.nc"

    if piSwitch:
        piCubeList = ants.load(inDir + piFiles)
        print()
        print("piCubeList:")
        print()
        print(piCubeList)
        print()
        piCube = piCubeList[0] + piCubeList[1] + piCubeList[2] + piCubeList[3]
        print()
        print("piCube BEFORE regrdding:")
        print()
        print(piCube)
        print()

        piCube = regrid(piCube)
        print(piCube)
        piCube.var_name = "n_dep"
        piCube.long_name = "NITROGEN DEPOSITION (kgN/m2/s)"
        piCube.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi(stash)
        piCube.attributes["grid_staggering"] = 6
        print(piCube.coord("time").units)
        piCube.coord("time").units = cf_units.Unit(
            "days since 1850-01-01 00:00:00", calendar="360_day"
        )
        piCube.coord("time").points = piCube.coord("time").points * 30.0
        print(piCube.coord("time"))
        piCube.coord("time").bounds = None
        #        ants.utils.guess_bounds(piCube.coord('time'), strict=True)  # SKL commented out
        ants.utils.coord.guess_bounds(piCube.coord("time"), strict=True)  # SKL added
        print(piCube)
        print(piCube.__class__.__name__)
        piCube.data = np.clip(piCube.data, 0.0, float("Inf"))
        ants.save(piCube, pi_ancil_name, saver="ancil")
        print("made " + pi_ancil_name)
    #
    if histSwitch:
        histCubeList = ants.load(inDir + histFiles)
        print(histCubeList)
        histCube = histCubeList[0] + histCubeList[1] + histCubeList[2] + histCubeList[3]
        print(histCube)
        print("Fails here:")
        histCube = regrid(histCube)
        print(histCube)
        histCube.var_name = "n_dep"
        histCube.long_name = "NITROGEN DEPOSITION (kgN/m2/s)"
        histCube.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi(stash)
        histCube.attributes["grid_staggering"] = 6
        print(histCube.coord("time").units)
        histCube.coord("time").units = cf_units.Unit(
            "days since 1850-01-01 00:00:00", calendar="360_day"
        )
        histCube.coord("time").points = histCube.coord("time").points * 30.0
        print(histCube.coord("time"))
        histCube.coord("time").bounds = None
        #       ants.utils.guess_bounds(histCube.coord('time'), strict=True) # SKL commented out
        ants.utils.coord.guess_bounds(histCube.coord("time"), strict=True)
        print(histCube)
        print(histCube.__class__.__name__)
        histCube.data = np.clip(histCube.data, 0.0, float("Inf"))
        ants.save(histCube, hist_ancil_name, saver="ancil")
        print("made " + hist_ancil_name)


def regrid(inCube):  # ,gridCube):
    """
    regrid the data onto the HG3 n96e grid
    new grid is specified by "gridFile"
    gridFile is a fractional area ancillary
    gridFile is not passed into the routine because the domain decomposition can't handle cubes with different shapes
    """
    gridCube = ants.load_cube(gridFile, iris.AttributeConstraint(STASH="m01s00i030"))
    inCube.coord("latitude").coord_system = gridCube.coord("latitude").coord_system
    inCube.coord("longitude").coord_system = gridCube.coord("longitude").coord_system
    inCube.coord("longitude").circular = True
    gridCube.coord("longitude").circular = True
    #    help(ants.utils.cube)

    #    ants.utils.guess_horizontal_bounds(inCube) # SKL commented out
    #    ants.utils.guess_horizontal_bounds(gridCube) # SKL commented out
    #    ants.utils.coord.guess_horizontal_bounds(inCube) # SKL added
    #    ants.utils.coord.guess_horizontal_bounds(gridCube) # SKL added
    ants.utils.cube.guess_horizontal_bounds(inCube)  # SKL added
    ants.utils.cube.guess_horizontal_bounds(gridCube)  # SKL added

    outCube = inCube.regrid(gridCube, iris.analysis.AreaWeighted())
    return outCube


if __name__ == "__main__":
    main()
