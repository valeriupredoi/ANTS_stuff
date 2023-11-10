"""
plot the original source files
against what UKCA has read in and processed: Test_RCP_Full.dat
the latter has got the conversion factors

T. Kuhlbrodt 23/11/17
"""

import glob
import os

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

SOURCE = "/group_workspaces/jasmin2/tids/CMIP6_ANCIL/data/inputs4MIPs_2017-06-08/UoM/GHGConcentrations/CMIP/UoM-CMIP-1-2-0/yr/"  # noqa
PART1 = "mole_fraction_of_"
PART2 = "_in_air"
PART3 = "gr1-GMNHSH/v20160902/"
ZERO = 0.0
DUM = -1.0e6
YEARDAY = 365.24

START = 1850
END = 2015
UKCAOUT = (
    "/group_workspaces/jasmin2/tids/users/till/RCP/UKCA_output/Test_RCP_Full_0line.dat"  # noqa
)

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
    ukcafull = np.genfromtxt(fname, dtype=float)
    ukcaout = ukcafull[::30, :]
    print("Shape:", np.shape(ukcaout))

    # this reads the header
    ukcaheader = []
    with open(fname, "r") as inp:
        for counter, line in enumerate(inp):
            if counter == 6:
                break
            ukcaheader.append(line.split())

    print("Shape:", np.shape(ukcaheader))

    for counter, value in enumerate(ukcaheader):
        print("Shape: ", np.size(value))
        print("line " + str(counter) + ":", value)

    species = ukcaheader[3][1:]
    convfac = np.asarray(ukcaheader[4][2:]).astype(float)
    # FLAKE8: unused variable
    # column = np.asarray(ukcaheader[5][2:]).astype(int)

    # now read the original source data fron netCDF
    start = START
    end = END
    # load files
    print("now the netCDF files!")
    molef = []
    units = []
    times = []
    for gas in GASES_NC:
        factor = 1.0
        varname = PART1 + STDIC[gas][0] + PART2
        print("VAR: ", varname)
        # FLAKE8: unused variable
        # path = os.path.join(SOURCE, varname, PART3, "*.nc")
        fname = glob.glob(os.path.join(SOURCE, varname, PART3, "*.nc"))
        if len(fname) != 1:
            raise Exception("Either too many or few input nc files " +
                            gas["name"])
        #   add some exceptions
        if gas[0:4] == "DUMM":
            factor = DUM
        if gas == "HALON1202":
            factor = ZERO

        fh = Dataset(fname[0], mode="r")
        molef.append(fh.variables[varname][start:end, 0] * factor)
        times.append(fh.variables["time"][start:end])
        units.append(fh.variables[varname].units)
        fh.close()

    size = np.shape(molef)
    print("number of variables ", size)

    taxis = times[0]

    # now do the plots
    # cycle through GASES
    plt.figure(figsize=(2 * 11.69, 2 * 8.27))
    print("PLOTTING NOW")
    plt.suptitle("netCDF in: green, UKCA out: red")
    for counter, spec in enumerate(GASES):

        ind = GASES_NC.index(spec)
        print("ind: ", ind)
        rcp = molef[:][ind]
        print("plot shape: ", np.shape(rcp[:]))
        unit = float(units[ind])
        print("unit: ", unit)

        ukcaind = species.index(MATCH[spec][0])
        print("ukcaind: ", ukcaind)
        print("UKCA species: ", species[ukcaind])
        conv = convfac[ukcaind - 1]
        print("convfac: ", conv)
        ukca = ukcaout[:350 * 12, ukcaind] / conv / unit
        #        print "UKCA time series: ", ukca

        plt.subplot(4, 5, counter + 1)
        plt.plot(taxis / YEARDAY, rcp[:], "g+", ukcaout[:350 * 12, 0], ukca,
                 "r-")
        plt.title(spec + " [" + units[ind] + "]")
    #        plt.xlabel("years")

    answer = input("Save this Figure? ")
    if answer == ("y" or "yes"):
        plt.savefig("./ghg_rcp_plot_allgases_nc_plus", dpi=350)
    else:
        plt.show()


if __name__ == "__main__":
    main()
