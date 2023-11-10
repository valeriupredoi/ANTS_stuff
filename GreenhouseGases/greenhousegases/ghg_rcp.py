import glob
import os

import numpy as np
from netCDF4 import Dataset


def __doc__():
    fstr = f"""
reformat the CMIP6 GHG data suitable for UKCA
use the RCP ASCII format and mole fractions
no conversion or interpolation needed

* read mole fractions, times, units
* put everything in a big list
* then write to ASCII format
  * three blanks between entries
  *  a blueprint for the data format is /group_workspaces/jasmin2/tids/users/till/RCP/trgas_rcp_scenario.dat

T. Kuhlbrodt 15/8/17

Below the variable names in the input4MIP data set
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

mole_fraction_of_c2f6_in_air                  mole_fraction_of_halon2402_in_air
mole_fraction_of_c3f8_in_air                  mole_fraction_of_hcfc141b_in_air
mole_fraction_of_c4f10_in_air                 mole_fraction_of_hcfc142b_in_air
mole_fraction_of_c5f12_in_air                 mole_fraction_of_hcfc22_in_air
mole_fraction_of_c6f14_in_air                 mole_fraction_of_hfc125_in_air
mole_fraction_of_c7f16_in_air                 mole_fraction_of_hfc134a_in_air
mole_fraction_of_c8f18_in_air                 mole_fraction_of_hfc134aeq_in_air
mole_fraction_of_c_c4f8_in_air                mole_fraction_of_hfc143a_in_air
mole_fraction_of_carbon_dioxide_in_air        mole_fraction_of_hfc152a_in_air
mole_fraction_of_carbon_tetrachloride_in_air  mole_fraction_of_hfc227ea_in_air
mole_fraction_of_cf4_in_air                   mole_fraction_of_hfc236fa_in_air
mole_fraction_of_cfc113_in_air                mole_fraction_of_hfc23_in_air
mole_fraction_of_cfc114_in_air                mole_fraction_of_hfc245fa_in_air
mole_fraction_of_cfc115_in_air                mole_fraction_of_hfc32_in_air
mole_fraction_of_cfc11_in_air                 mole_fraction_of_hfc365mfc_in_air
mole_fraction_of_cfc11eq_in_air               mole_fraction_of_hfc4310mee_in_air
mole_fraction_of_cfc12_in_air                 mole_fraction_of_methane_in_air
mole_fraction_of_cfc12eq_in_air               mole_fraction_of_methyl_bromide_in_air
mole_fraction_of_ch2cl2_in_air                mole_fraction_of_methyl_chloride_in_air
mole_fraction_of_ch3ccl3_in_air               mole_fraction_of_nf3_in_air
mole_fraction_of_chcl3_in_air                 mole_fraction_of_nitrous_oxide_in_air
mole_fraction_of_halon1211_in_air             mole_fraction_of_sf6_in_air
mole_fraction_of_halon1301_in_air             mole_fraction_of_so2f2_in_air

"""  # noqa

    return fstr


SOURCE = "/group_workspaces/jasmin2/tids/CMIP6_ANCIL/data/inputs4MIPs_2017-06-08/UoM/GHGConcentrations/CMIP/UoM-CMIP-1-2-0/yr/"  # noqa
PART1 = "mole_fraction_of_"
PART2 = "_in_air"
PART3 = "gr1-GMNHSH/v20160902/"
ZERO = 0.0
DUM = -1.0e6

START = 2014
END = 2015

TARGET = "./trgas_rcp_historical_" + str(START) + "_" + str(END) + ".dat"

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

OMAG = {"1.e-6": ["ppm"], "1.e-9": ["ppb"], "1.e-12": ["ppt"]}

GASES = [
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
    start = START
    end = END
    # load files
    molef = []
    units = []
    times = []
    for gas in GASES:
        factor = 1.0
        varname = PART1 + STDIC[gas][0] + PART2
        print("VAR: ", varname)
        # FLAKE8: unused variable
        # path = os.path.join(SOURCE, varname, PART3, "*.nc")
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
        fh.close()

    size = np.shape(molef)
    print("number of variables ", size[0])

    taxis = times[0]

    #   now write data to ASCII file
    fname = TARGET
    with open(fname, "w") as outp:
        #   the header
        outp.write("! historical GHG data for CMIP6 \n! source: " + SOURCE + " \n")
        outp.write(" &THISFILE_SPECIFICATIONS \n")
        outp.write(" THISFILE_DATACOLUMNS    =          " + str(size[0]) + ", \n")
        outp.write(" THISFILE_FIRSTYEAR      =        " + str(START) + ", \n")
        outp.write(" THISFILE_LASTYEAR       =        " + str(END - 1) + ", \n")
        outp.write(" THISFILE_ANNUALSTEPS    =           1, \n")
        outp.write(" THISFILE_FIRSTDATAROW   =          15, \n")
        outp.write(" THISFILE_UNITS          =  SEE ROW 13, \n")
        outp.write(" THISFILE_DATTYPE        =  RCPDAT \n/ \n")
        #   the data
        outp.write("   COLUMN:   ")
        for i in np.arange(size[0]):
            outp.write(
                "   %i" % (i + 1),
            )
        outp.write("\n   UNITS:   ")
        for counter, value in enumerate(units):
            outp.write(
                "   %s" % OMAG[value][0],
            )
        outp.write("\n   YEARS   ")
        for i in np.arange(size[0]):
            outp.write(
                "   %s" % GASES[i],
            )
        for j in np.arange(size[1]):
            outp.write(
                "\n   %i" % int(taxis[j] / 365.24),
            )
            for i in np.arange(size[0]):
                outp.write("    %7.4e" % molef[i][j])
        outp.write("\n")


if __name__ == "__main__":
    main()
