"""
reformat the CMIP6 GHG data suitable for the rose suite configuration file
convert from mole fraction into mass mixing ratio

things to do:

/ probably the metadata for the suite could be modified to have tick boxes to switch the hist. forcing on or off 
/ use the conversion factors from ukca_constants.F90 (15-7-16)
/ the time series are linearly interpolated because the UM expects 1st Jan values for each year. 
The input time series is global means from 1848 to 2014. The
output time series is interpolated to run from 1st Jan 1849 till 1st Jan 2015

T. Kuhlbrodt 2/12/16

"""

import glob
import os

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

# conversion factors from ukca_constants.F90
CO2_M = 1.5188
CH4_M = 0.5523
N2O_M = 1.5188
CFC_M = 4.1783
HFC_M = 3.5219

# Units from the netcdf files in the Source
CO2_U = 1.0e-6
CH4_U = 1.0e-9
N2O_U = 1.0e-9
CFC_U = 1.0e-12
HFC_U = 1.0e-12

CO2 = {
    "name": "CO2",
    "converfac": CO2_M,
    "units": CO2_U,
    "varname": "mole_fraction_of_carbon_dioxide_in_air",
}
CH4 = {
    "name": "CH4",
    "converfac": CH4_M,
    "units": CH4_U,
    "varname": "mole_fraction_of_methane_in_air",
}
N2O = {
    "name": "N2O",
    "converfac": N2O_M,
    "units": N2O_U,
    "varname": "mole_fraction_of_nitrous_oxide_in_air",
}
CFC12 = {
    "name": "CFC12",
    "converfac": CFC_M,
    "units": CFC_U,
    "varname": "mole_fraction_of_cfc12eq_in_air",
}
HFC = {
    "name": "HFC134A",
    "converfac": HFC_M,
    "units": HFC_U,
    "varname": "mole_fraction_of_hfc134aeq_in_air",
}
CFC11 = {"name": "CFC11", "converfac": "", "units": "", "varname": ""}
CFC113 = {"name": "CFC113", "converfac": "", "units": "", "varname": ""}
CFC114 = {"name": "CFC114", "converfac": "", "units": "", "varname": ""}
HCFC22 = {"name": "HCFC22", "converfac": "", "units": "", "varname": ""}
HFC125 = {"name": "HFC125", "converfac": "", "units": "", "varname": ""}
SO4 = {"name": "SO4", "converfac": "", "units": "", "varname": ""}

GASES = [CO2, CH4, N2O, CFC12, HFC]
GASES_NON_UM = [CFC11, CFC113, CFC114, HCFC22, HFC125, SO4]

SOURCEDIR = "/group_workspaces/jasmin2/tids/users/till/Source/03Oct/"


def interpol(gas):
    # linear interpolation
    gas_ipl = np.zeros(len(gas))
    gas_ipl[:-1] = 0.5 * (gas[1:] + gas[:-1])
    # linear extrapolation to 1st Jan 2015
    gas_ipl[-1] = gas[-1] + 0.5 * (gas[-1] - gas[-2])

    return gas_ipl


def main(gas_in):
    fname = glob.glob(os.path.join(SOURCEDIR, gas_in["name"], "*.nc"))
    print("FNAME", fname)
    if len(fname) != 1:
        raise Exception("Either too many or few input nc files " + gas["name"])
    fh = Dataset(fname[0], mode="r")
    # read only the first column (sector "GM", as opposed to the individual hemispheres)
    # read the last 167 years, i.e. from 1848
    gas = fh.variables[gas_in["varname"]][-167:, 0]
    fh.close()
    # convert volume mix ratio (i.e. mole fraction) to mass mixing ratio
    # use the units provided in the Source
    leng = len(gas)
    gas = gas * gas_in["converfac"] * gas_in["units"]
    # linear interpolation in time, plus an extrapolation to 1st Jan 2015
    gas = interpol(gas)
    # add 1849 to time axis
    time = np.arange(leng) + 1849

    # write out the time series in the rose suite configuration file format (csv with "=" and newlines thrown in)
    # fname=SOURCEDIR+GAS+"/"+GAS+".cnf"
    fname = os.path.join(SOURCEDIR, gas_in["name"], gas_in["name"] + ".cnf")
    with open(fname, "w") as outp:
        for line in range(leng / 6 + 1):
            outp.write("    =")
            for i in range(6):
                try:
                    outp.write("%7.4e," % gas[line * 6 + i])
                except:
                    pass
            outp.write("\n")

        for line in range(leng / 6 + 1):
            outp.write("    =")
            for i in range(6):
                try:
                    outp.write("%d," % time[line * 6 + i])
                except:
                    pass
            outp.write("\n")

    return gas, time


def write_rose_conf(gas_mmr):
    # now write an opt file for rose
    fname = os.path.join(SOURCEDIR, "global_ghg_cmip6.conf")
    print(gas_mmr)
    ntimes = len(gas_mmr[(GASES[0]["name"], "year")])
    print(ntimes)
    # write concentrations
    with open(fname, "w") as outp:
        outp.write("[namelist:clmchfcg]")
        outp.write("\n")
        for gas in GASES_NON_UM:
            outp.write("clim_fcg_levls_" + gas["name"].lower() + "=" +
                       str(ntimes) + "*-32768.0")
            outp.write("\n")
        for gas in GASES:
            for line in range(ntimes / 6 + 1):
                print(line)
                if line == 0:
                    outp.write("clim_fcg_levls_" + gas["name"].lower() + "=")
                else:
                    outp.write("    =")
                for i in range(6):
                    try:
                        outp.write("%7.4e," %
                                   gas_mmr[(gas["name"], "mmr")][line * 6 + i])
                    except:
                        pass
                outp.write("\n")
        # write no years
        for gas in GASES_NON_UM:
            outp.write("clim_fcg_nyears_" + gas["name"].lower() + "=0")
            outp.write("\n")
        for gas in GASES:
            outp.write("clim_fcg_nyears_" + gas["name"].lower() + "=" +
                       str(ntimes))
            outp.write("\n")
        # write rates
        for gas in GASES_NON_UM:
            outp.write("clim_fcg_rates_" + gas["name"].lower() + "=" +
                       str(ntimes) + "*-32768.0")
            outp.write("\n")
        for gas in GASES:
            outp.write("clim_fcg_rates_" + gas["name"].lower() + "=" +
                       str(ntimes) + "*-32768.0")
            outp.write("\n")
        # write years
        for gas in GASES_NON_UM:
            outp.write("clim_fcg_years_" + gas["name"].lower() + "=" +
                       str(ntimes) + "*-32768")
            outp.write("\n")
        loop_step = 12
        for gas in GASES:
            for line in range(ntimes / loop_step + 1):
                print(line)
                if line == 0:
                    outp.write("clim_fcg_years_" + gas["name"].lower() + "=")
                else:
                    outp.write("    =")
                for i in range(loop_step):
                    try:
                        outp.write("%d," %
                                   gas_mmr[(gas["name"],
                                            "year")][line * loop_step + i])
                    except:
                        pass
                outp.write("\n")


if __name__ == "__main__":
    gas_mmr = {}
    for gas in GASES:
        mmr, year = main(gas)
        gas_mmr[(gas["name"], "mmr")] = mmr
        gas_mmr[(gas["name"], "year")] = year

    write_rose_conf(gas_mmr)
