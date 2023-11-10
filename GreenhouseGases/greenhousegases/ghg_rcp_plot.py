"""
plot the output of ghg_rcp.py: trgas_rcp_hitorical.dat
against what UKCA has read in and processed: Test_RCP_Full.dat
the latter has got the conversion factors

T. Kuhlbrodt 21/9/17

trgas_rcp_historical.dat has got these variables:
  [ YEARS]      DUMM0   CO2   CH4   N2O   HFC125   HFC134a   CFC_11   CFC_12   DUMM1   CFC_113   CFC_114   CFC_115   CARB_TET   MCF   HCFC_22   HCFC_141B   HCFC_142B   HALON1211   HALON1202   HALON1301   HALON2402   CH3BR   CH3CL   DUMM2

Test_RCP_Full.dat has got these variables: 
#      SPECIES                CFCl3               CF2Cl2           CF2ClCFCl2           CF2ClCF2Cl             CF2ClCF3                 CCl4               MeCCl3               CHF2Cl              MeCF2Cl              MeCFCl2              CF2ClBr               CF2Br2                CF3Br           CF2BrCF2Br                 MeBr                 MeCl                  N2O                  CH4                   H2                   N2                  CO2               CH2Br2                CHBr3              CF3CHF2              CH2FCF3                  COS

In this table I am matching the variables from trgas to those from Test_RCP_Full:

trgas | Test_RCP_Full

CO2   | CO2
CH4   | CH4
N2O   | N2O
HFC125 | CF3CHF2
HFC134a | CH2FCF3
CFC_11 |  CFCl3
CFC_12 | CF2Cl2
CFC_113 | CF2ClCFCl2
CFC_114 | CF2ClCF2Cl
CFC_115 | CF2ClCF3
CARB_TET | CCl4 
MCF | MeCCl3
HCFC_22 | CHF2Cl
HCFC_141B | MeCFCl2
HCFC_142B | MeCF2Cl
HALON1211 |  CF2ClBr
HALON1202 | CF2Br2
HALON1301 | CF3Br
HALON2402 | CF2BrCF2Br
CH3BR | MeBr
CH3CL | MeCl

"""

import glob
import os

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

SOURCE = "./trgas_rcp_historical.dat"
# UKCAOUT='/group_workspaces/jasmin2/tids/users/till/RCP/UKCA_output/Test_RCP_Short.dat'
UKCAOUT = (
    "/group_workspaces/jasmin2/tids/users/till/RCP/UKCA_output/Test_RCP_Full_0line.dat"
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

UNITS = {"ppm": [1.0e-6], "ppb": [1.0e-9], "ppt": [1.0e-12]}

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
    column = np.asarray(ukcaheader[5][2:]).astype(int)

    # now read the RCP input file
    fname = SOURCE
    # this reads the actual data
    rcpin = np.genfromtxt(fname, dtype=float, skip_header=14)
    print("RCPin shape: ", np.shape(rcpin))

    # this reads the header
    rcpheader = []
    with open(fname, "r") as inp:
        for counter, line in enumerate(inp):
            if counter == 16:
                break
            rcpheader.append(line.split())

    print("Shape:", np.shape(rcpheader))

    for counter, value in enumerate(rcpheader):
        print("Shape: ", np.size(value))
        print("line " + str(counter) + ":", value)

    rcpspecies = rcpheader[13][:]
    rcpcolumn = np.asarray(rcpheader[11][1:]).astype(int)
    rcpunits = rcpheader[12][:]

    # now do the plots
    # cycle through GASES
    plt.figure(figsize=(2 * 11.69, 2 * 8.27))
    plt.suptitle("RCP in: blue, UKCA out: red")
    for counter, spec in enumerate(GASES):

        ind = rcpspecies.index(spec)
        print("ind: ", ind)
        rcp = rcpin[:, ind]
        print("plot shape: ", np.shape(rcp[:]))
        unitstr = rcpunits[ind]
        unit = UNITS[unitstr]
        print("unit: ", unit)

        ukcaind = species.index(MATCH[spec][0])
        print("ukcaind: ", ukcaind)
        print("UKCA species: ", species[ukcaind])
        conv = convfac[ukcaind - 1]
        print("convfac: ", conv)
        ukca = ukcaout[: 350 * 12, ukcaind] / conv / unit
        #        print "UKCA time series: ", ukca

        plt.subplot(4, 5, counter + 1)
        plt.plot(rcpin[:, 0], rcp[:], "b", ukcaout[: 350 * 12, 0], ukca, "r-")
        plt.title(spec + " [" + unitstr + "]")
    #        plt.xlabel("years")

    answer = input("Save this Figure? ")
    if answer == ("y" or "yes"):
        plt.savefig("./ghg_rcp_plot_allgases", dpi=350)
    else:
        plt.show()


if __name__ == "__main__":
    main()
