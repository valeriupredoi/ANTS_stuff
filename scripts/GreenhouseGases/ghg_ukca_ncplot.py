#!/usr/bin/env python
"""
plot the original source files
against what GHG_UKCA has read in and processed: trgas_ssp370_2015_2200.dat 

T. Kuhlbrodt 29/11/18
"""

import glob
import os

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

SOURCE = "/gws/nopw/j04/cmip6_prep_vol1/cmip6_ancils/data/inputs4MIPs_2018-11-01/CMIP6/ScenarioMIP/UoM/UoM-AIM-ssp370-1-2-0/atmos/yr/"
PART1 = "mole_fraction_of_"
PART2 = "_in_air"
PART3 = "gr1-GMNHSH/v20180611"
ZERO = 0.0
DUM = -1.0e6
YEARDAY = 365.0

START = 2015
END = 2200
UKCAOUT = "/gws/nopw/j04/cmip6_prep_vol1/cmip6_ancils/users/till/Out/trgas_ssp370_2015_2200.dat"

MATCH = {
    "CO2": ["CO2"],
    "CH4": ["CH4"],
    "N2O": ["N2O"],
    "HFC125": ["CF3CHF2"],
    "HFC134a": ["CH2FCF3"],
    "CFC_11": ["CFCl3"],
    "CFC_12": ["CF2Cl2"],
    "CFC_113": ["CF2ClCFCl2"],
    "CFC_114": ["CF2ClCF2Cl"],
    "CFC_115": ["CF2ClCF3"],
    "CARB_TET": ["CCl4"],
    "MCF": ["MeCCl3"],
    "HCFC_22": ["CHF2Cl"],
    "HCFC_141B": ["MeCFCl2"],
    "HCFC_142B": ["MeCF2Cl"],
    "HALON1211": ["CF2ClBr"],
    "HALON1202": ["CF2Br2"],
    "HALON1301": ["CF3Br"],
    "HALON2402": ["CF2BrCF2Br"],
    "CH3BR": ["MeBr"],
    "CH3CL": ["MeCl"],
}

STDIC = {
    "CO2": ["carbon_dioxide"],
    "CH4": ["methane"],
    "N2O": ["nitrous_oxide"],
    "HFC134a": ["hfc134a"],
    "CFC_12": ["cfc12"],
    "DUMM0": ["cfc12"],
    "HFC125": ["hfc125"],
    "CFC_11": ["cfc11"],
    "DUMM1": ["cfc12"],
    "CFC_113": ["cfc113"],
    "CFC_114": ["cfc114"],
    "CFC_115": ["cfc115"],
    "CARB_TET": ["carbon_tetrachloride"],
    "MCF": ["ch3ccl3"],
    "HCFC_22": ["hcfc22"],
    "HCFC_141B": ["hcfc141b"],
    "HCFC_142B": ["hcfc142b"],
    "HALON1211": ["halon1211"],
    "HALON1202": ["halon1211"],
    "HALON1301": ["halon1301"],
    "HALON2402": ["halon2402"],
    "CH3BR": ["methyl_bromide"],
    "CH3CL": ["methyl_chloride"],
    "DUMM2": ["cfc12"],
}

GASES = [
    "CO2",
    "CH4",
    "N2O",
    "HFC125",
    "HFC134a",
    "CFC_11",
    "CFC_12",
    "CFC_113",
    "CFC_114",
    "CFC_115",
    "CARB_TET",
    "MCF",
    "HCFC_22",
    "HCFC_141B",
    "HCFC_142B",
    "HALON1211",
    "HALON1301",
    "HALON2402",
    "CH3BR",
    "CH3CL",
]

GASES_NC = [
    "DUMM0",
    "CO2",
    "CH4",
    "N2O",
    "HFC125",
    "HFC134a",
    "CFC_11",
    "CFC_12",
    "DUMM1",
    "CFC_113",
    "CFC_114",
    "CFC_115",
    "CARB_TET",
    "MCF",
    "HCFC_22",
    "HCFC_141B",
    "HCFC_142B",
    "HALON1211",
    "HALON1202",
    "HALON1301",
    "HALON2402",
    "CH3BR",
    "CH3CL",
    "DUMM2",
]


def main():
    # first read UKCAOUT
    # line 4 is species
    # line 5 is conversion factor
    # from line 7: years

    fname = UKCAOUT
    # this reads the actual data
    ukcafull = np.genfromtxt(fname, dtype=float, skip_header=14)
    ukcaout = ukcafull[:, :]
    print("Shape:", np.shape(ukcaout))

    # this reads the header
    ukcaheader = []
    with open(fname, "r") as inp:
        for counter, line in enumerate(inp):
            if counter == 14:
                break
            ukcaheader.append(line.split())

    print("Shape:", np.shape(ukcaheader))

    for counter, value in enumerate(ukcaheader):
        print("Shape: ", np.size(value))
        print("line " + str(counter) + ":", value)

    species = ukcaheader[13][1:]
    print("Species: ", species)
    #    convfac=np.asarray(ukcaheader[4][2:]).astype(float)
    column = np.asarray(ukcaheader[11][1:]).astype(int)

    # now read the original source data from netCDF
    start = START - 2015
    end = END - 2015
    # load files
    print("now the netCDF files!")
    molef = []
    units = []
    times = []
    gasspec = []
    for gas in GASES_NC:
        factor = 1.0
        varname = PART1 + STDIC[gas][0] + PART2
        print("VAR: ", varname)
        path = os.path.join(SOURCE, varname, PART3, "*.nc")
        fname = glob.glob(os.path.join(SOURCE, varname, PART3, "*.nc"))
        if len(fname) != 1:
            raise Exception("Either too many or few input nc files " + gas["name"])
        #   add some exceptions
        if gas[0:4] == "DUMM":
            factor = DUM
        if gas == "HALON1202":
            factor = ZERO

        fh = Dataset(fname[0], mode="r")
        molef.append(fh.variables[varname][start:end, 0] * factor)
        times.append(fh.variables["time"][start:end])
        units.append(fh.variables[varname].units)
        gasspec.append(gas)
        fh.close()

    size = np.shape(molef)
    print("number of variables ", size)

    taxis = times[0]

    # now do the plots
    # cycle through GASES
    plt.figure(figsize=(2 * 11.69, 2 * 8.27))
    print("PLOTTING NOW")
    plt.suptitle("netCDF in: green, UKCA ancil: red \n Scenario: UoM-AIM-ssp370-1-2-0")
    for counter, spec in enumerate(GASES):

        ind = GASES_NC.index(spec)
        print("ind: ", ind)
        print("now plotting gas: ", spec, " i.e. ", gasspec[
            ind
        ], " from nc file and gas ", GASES_NC[ind], "from UKCA ASCII file")
        rcp = molef[:][ind]
        #        print "plot shape: ", np.shape(rcp[:])
        # units in source files are wrong
        if spec == "CO2":
            unit = 1.0e-6
        elif spec == "N2O" or spec == "CH4":
            unit = 1.0e-9
        else:
            unit = 1.0e-12
        print("unit: ", unit)

        ukcaind = species.index(spec)
        print("ukcaind: ", ukcaind)
        print("UKCA species: ", species[ukcaind])

        plt.subplot(4, 5, counter + 1)
        plt.plot(
            taxis / YEARDAY + 1850,
            rcp[:],
            "g+",
            ukcaout[:, 0],
            ukcaout[:, ukcaind + 1],
            "r-",
        )
        plt.title(spec + " [" + str(unit) + "]")
    #        plt.xlabel("years")

    answer = input("Save this Figure? ")
    if answer == "y" or answer == "yes":
        plt.savefig("./ghg_ukca_plot_allgases_nc_plus", dpi=350)
    else:
        plt.show()


if __name__ == "__main__":
    main()
