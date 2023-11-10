#!/usr/bin/env python
import argparse
import glob
import os
from warnings import warn

import numpy as np
from netCDF4 import Dataset


def __doc__():
    fstr = f"""
        reformat the CMIP6 GHG data suitable for the rose suite configuration file
        convert from mole fraction into mass mixing ratio
        a version with a command line interface

        interpolate the source data in time
        if this fails for specific years (around 2015AD for instance), then:
        * at the beginning, pad with the first available year from the source
        * at the end, extrapolate from the source 

        Location of HISTORICAL source data: 
        /gws/nopw/j04/cmip6_prep_vol1/cmip6_ancils/data/inputs4MIPs_2018-11-01/CMIP6/CMIP/UoM/UoM-CMIP-1-2-0/atmos/yr/mole-fraction-of-carbon-dioxide-in-air/gr1-GMNHSH/v20160830/

        Location of SCENARIOMIP source data:
        /gws/nopw/j04/cmip6_prep_vol1/cmip6_ancils/data/inputs4MIPs_2018-11-01/CMIP6/ScenarioMIP/UoM/UoM-AIM-ssp370-1-2-0/atmos/yr/mole_fraction_of_carbon_dioxide_in_air/gr1-GMNHSH/v20180611

        I'm defining the parameters for the in-line command:
        $> source='/gws/nopw/j04/cmip6_prep_vol1/cmip6_ancils/data/inputs4MIPs_2018-11-01/CMIP6/ScenarioMIP/UoM/UoM-IMAGE-ssp126-1-2-0/atmos/yr'
        $> target='./Outnow/glo_ghg_SSP126_'

        and then I can call the routine like so:
        $> ./GHG_radiation.py -s $source -o $target -p "scenarioMIP" -b 2015 -e 2101

        For the SSP585 data: source='/gws/nopw/j04/cmip6_prep_vol1/cmip6_ancils/data/inputs4MIPs_2018-11-01/CMIP6/ScenarioMIP/UoM/UoM-REMIND-MAGPIE-ssp585-1-2-0/atmos/yr'
        $> target='./Outnow/glo_ghg_SSP585_'
        Alternatively, for the historical data:
        $> source='/gws/nopw/j04/cmip6_prep_vol1/cmip6_ancils/data/inputs4MIPs_2018-11-01/CMIP6/CMIP/UoM/UoM-CMIP-1-2-0/atmos/yr'
        $> target='./Outnow/glo_ghg_hist_'

        and then I can call the routine like so:
        $> ./GHG_scenario.py -s $source -o $target -p "historical" -b 1850 -e 2015

        T. Kuhlbrodt 14/2/19
        """  # noqa

    return fstr


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
    "varname": "carbon_dioxide_GM",
}
CH4 = {
    "name": "CH4",
    "converfac": CH4_M,
    "units": CH4_U,
    "varname": "methane_GM"
}
N2O = {
    "name": "N2O",
    "converfac": N2O_M,
    "units": N2O_U,
    "varname": "nitrous_oxide_GM"
}
CFC12 = {
    "name": "CFC12",
    "converfac": CFC_M,
    "units": CFC_U,
    "varname": "cfc12eq_GM"
}
HFC = {
    "name": "HFC134A",
    "converfac": HFC_M,
    "units": HFC_U,
    "varname": "hfc134aeq_GM"
}
CFC11 = {"name": "CFC11", "converfac": "", "units": "", "varname": ""}
CFC113 = {"name": "CFC113", "converfac": "", "units": "", "varname": ""}
CFC114 = {"name": "CFC114", "converfac": "", "units": "", "varname": ""}
HCFC22 = {"name": "HCFC22", "converfac": "", "units": "", "varname": ""}
HFC125 = {"name": "HFC125", "converfac": "", "units": "", "varname": ""}
SO4 = {"name": "SO4", "converfac": "", "units": "", "varname": ""}

GASES = [CO2, CH4, N2O, CFC12, HFC]
GASES_NON_UM = [CFC11, CFC113, CFC114, HCFC22, HFC125, SO4]

# boundary between historical and scenarioMIP simulations
WATERSHED = 2015


def interpolate(array):
    # linear interpolation in time. Returned array has a length reduced by one
    return 0.5 * (array[1:] + array[:-1])


def ghg_ssp(gas_in, start, end, source_files, project, data_version):
    # nc file calendar parameters
    if project == "historical":
        FIRST = 0
        DAYS = 365.24
        if data_version is None:
            data_version = "v20160830"
    elif project == "scenarioMIP":
        FIRST = 1850
        DAYS = 365.0
        if data_version is None:
            data_version = "v20180611"
    else:
        print("MIP project needs to be either 'historical' or 'scenarioMIP'.")
        exit()
    # sanity check: historical end year can't be greater than 2015;
    # likewise scenarioMIP start year can't be smaller than 2015
    if project == "historical":
        if end > WATERSHED:
            warn((
                "'end' argument {} is later than final year "
                "in source file. Resetting to {}"
            ).format(end, WATERSHED))
            end = WATERSHED
    elif project == "scenarioMIP":
        if start < WATERSHED:
            warn((
                "'start' argument {} is earlier than start year "
                "in source file. Resetting to {}"
            ).format(start, WATERSHED))
            start = WATERSHED
    # fixed parts of the source filename
    PART1 = "mole_fraction_of_"
    PART2 = "_in_air"
    PART3 = os.path.join("gr1-GMNHSH", data_version)

    variname = PART1 + gas_in["varname"][:-3] + PART2
    if project == "historical":  # turn underscores into dashes
        variname_mod = PART1 + gas_in["varname"][:-3] + PART2
        variname_mod = variname_mod.replace("_", "-")
        fname = glob.glob(
            os.path.join(source_files, variname_mod, PART3, "*.nc"))
    elif project == "scenarioMIP":
        fname = glob.glob(os.path.join(source_files, variname, PART3, "*.nc"))
    print((fname[0]), end=" ")
    if len(fname) != 1:
        raise Exception("Either too many or no input nc files " +
                        gas_in["name"])

    fh = Dataset(fname[0], mode="r")
    # read in the data
    gas = fh.variables[variname][:, 0]
    time = fh.variables["time"][:] / DAYS + FIRST
    fh.close

    # convert volume mix ratio (i.e. mole fraction) to mass mixing ratio
    # use the units provided in the Source
    # leng = len(gas)  # FLAKE8: unused variable
    gas = gas * gas_in["converfac"] * gas_in["units"]

    if end:
        # throw an error if start >= end
        if start >= end:
            print("Starting year needs to be earlier than finishing year.")
            print("You've got: start = ", start, "and end = ", end)
            exit()
        # linear interpolation in time
        gas_intp = interpolate(gas)
        time_intp = interpolate(time)
        # set the actual start and end year to be cut out
        # don't add an extra year at the end - if extra years
        # are needed the user has to request it
        start_read = start
        end_read = end

        # check if year (start) is available in the interpolated time series
        if time_intp[0] > start:
            start_read = start_read + 1
            # convert time array to int, this is what's needed
            # for the rose-app.conf format
            time_intp = time_intp.astype(int)
            chunk = np.where((time_intp >= start_read)
                             & (time_intp <= end_read))
            gas_out = gas_intp[chunk]
            time_out = time_intp[chunk]
            # padding at the beginning if year (start - 1)
            # is not available in the source.
            # do this from the source data (not from the interpolated data)
            # FLAKE8: unused variable
            # chunk_extrap = np.where(np.isclose(time, start + 0.5))
            gas_extrap = gas[chunk]
            gas_out = np.insert(gas_out, 0, gas_extrap[0])
            time_out = np.insert(time_out, 0, time_out[0] - 1)

        # check if year (end) is available in the interpolated time series
        if time_intp[-1] < end:
            end_read = end - 1
            # convert time array to int
            # this is what's needed for the rose-app.conf format
            time_intp = time_intp.astype(int)
            chunk = np.where((time_intp >= start_read)
                             & (time_intp <= end_read))
            gas_out = gas_intp[chunk]
            time_out = time_intp[chunk]
            # extrapolation at the end, from the source data
            # (not from the interpolated data)
            # chunk_extrap = np.where((time <= end - 0.49)
            #                        & (time >= end - 1.51))
            # FLAKE8: unused variable
            gas_extrap = gas[chunk]
            gas_out = np.insert(
                gas_out, -1,
                gas_extrap[-1] + 0.5 * (gas_extrap[-1] - gas_extrap[-2]))
            time_out = np.insert(time_out, -1, time_out[-1] + 1)

    # now the case where end==None, that is, just a single year
    else:
        # convert time array to int
        # is what's needed for the rose-app.conf format
        time = time.astype(int)
        chunk = np.where(time == start)
        gas_out = gas[chunk]
        time_out = time[chunk]

    return gas_out, time_out


def write_rose_conf(gas_mmr, start, end, output_file, project):
    #   now write an opt file for rose
    if end is None:
        output = output_file + str(start) + ".conf"
    else:
        output = output_file + str(start) + "_" + str(end) + ".conf"
    #   print gas_mmr
    ntimes = len(gas_mmr[(GASES[0]["name"], "year")])
    # write concentrations
    with open(output, "w") as outp:
        outp.write("[namelist:clmchfcg]")
        outp.write("\n")
        for gas in GASES_NON_UM:
            outp.write("clim_fcg_levls_" + gas["name"].lower() + "=" +
                       str(ntimes) + "*-32768.0")
            outp.write("\n")
        loop_step = 6
        for gas in GASES:
            for line in range(ntimes / loop_step + 1):
                #                print line
                if line == 0:
                    outp.write("clim_fcg_levls_" + gas["name"].lower() + "=")
                else:
                    outp.write("    =")
                for i in range(loop_step):
                    count = line * loop_step + i
                    #  last entry in list shouldn't be followed by a comma
                    # (makes rose GUI throw an error)
                    if count < (ntimes - 1):
                        outp.write("%7.4e," %
                                   gas_mmr[(gas["name"], "mmr")][count])
                    else:
                        try:
                            outp.write("%7.4e" %
                                       gas_mmr[(gas["name"], "mmr")][count])
                        except RuntimeError as exc:
                            print(exc)
                            break
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
                #                print line
                if line == 0:
                    outp.write("clim_fcg_years_" + gas["name"].lower() + "=")
                else:
                    outp.write("    =")
                for i in range(loop_step):
                    count = line * loop_step + i
                    if count < (ntimes - 1):
                        outp.write("%d," %
                                   gas_mmr[(gas["name"], "year")][count])
                    else:
                        try:
                            outp.write("%d" %
                                       gas_mmr[(gas["name"], "year")][count])
                        except RuntimeError as exc:
                            print(exc)
                            break
                outp.write("\n")


def _process(start, end, source_files, output_file, project, data_version):
    # Distinguish between timeslice/timeseries based on whether an end
    # argument is provided:
    gas_mmr = {}
    #    print "GASES: ", GASES
    for gas in GASES:
        mmr, year = ghg_ssp(gas, start, end, source_files[0], project,
                            data_version)
        gas_mmr[(gas["name"], "mmr")] = mmr
        gas_mmr[(gas["name"], "year")] = year

    write_rose_conf(gas_mmr, start, end, output_file, project)


if __name__ == "__main__":
    # __doc__ is the module docstring.
    arg_parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    args = arg_parser.parse_args()
    _process(args.begin, args.end, args.sources, args.output, args.project,
             args.data_version)


def __example__():
    fstr = f"""
        Example usages:

        * The help interface:

        $ ./example.py -h

        usage: example.py [-h] [-e END] [-n NAME] -b BEGIN -o OUTPUT -s SOURCES [SOURCES ...] [-f FOO]

        This module docstring will be passed to the ArgumentParser to provide usage
        instructions for the application.  Notice how the help text is populated with
        little effort.

        optional arguments:
          -h, --help            show this help message and exit
          -e END, --end END     Optional end year for the processing. If this is used, output is a
                                timeseries from the 'begin' year to this year, inclusive of both the
                                'begin' year and the 'end' year. Otherwise, the output will be a
                                timeslice in the 'begin' year only.
          -n NAME, --name NAME  Name of the species to be preprocessed.
          -f FOO, --foo FOO     An arbitrary new argument.

        required arguments:
          -b BEGIN, --begin BEGIN
                                Start year for the processing.
          -o OUTPUT, --output OUTPUT
                                Filename to write the result to.
          -s SOURCES [SOURCES ...], --sources SOURCES [SOURCES ...]
                                Filename(s) to read the source data from.

        * Generating a timeslice from two source files:

        $ ./example.py -b 1850 -n SO2 -s foo bar -o baz
        Process a timeslice for species SO2 from ['foo', 'bar'] files for year 1850, writing result to baz.

        * Generating a timeseries:

        $ ./example.py -b 1850 -e 2014 -n SO2 -s foo bar -o baz
        Process a timeseries for species SO2 from ['foo', 'bar'] files beween 1850 and 2014, writing result to baz.

        """  # noqa

    return fstr
