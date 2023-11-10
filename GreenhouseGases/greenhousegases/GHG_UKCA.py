#!/usr/bin/env python
import argparse
import glob
import os
from warnings import warn

import numpy as np
from netCDF4 import Dataset


def __doc__():
    fstr = f"""
        reformat the CMIP6 GHG data suitable for UKCA
        use the RCP ASCII format and mole fractions
        no conversion or interpolation needed

        / read mole fractions, times, units
        / put everything in a big list
        / then write to ASCII format
          * three blanks between entries
          *  a blueprint for the data format is /group_workspaces/jasmin2/tids/users/till/RCP/trgas_rcp_scenario.dat
  
       * use with command line arguments like GHG_SSP.py
       * new input (Source) data needed

       Location of HISTORICAL source data: 
       /gws/nopw/j04/cmip6_prep_vol1/cmip6_ancils/data/inputs4MIPs_2018-11-01/CMIP6/CMIP/UoM/UoM-CMIP-1-2-0/atmos/yr/mole-fraction-of-carbon-dioxide-in-air/gr1-GMNHSH/v20160830/

       Location of SCENARIOMIP source data:
       /gws/nopw/j04/cmip6_prep_vol1/cmip6_ancils/data/inputs4MIPs_2018-11-01/CMIP6/ScenarioMIP/UoM/UoM-AIM-ssp370-1-2-0/atmos/yr/mole_fraction_of_carbon_dioxide_in_air/gr1-GMNHSH/v20180611

      I'm defining the parameters for the in-line command:
      $> source='/gws/nopw/j04/cmip6_prep_vol1/cmip6_ancils/data/inputs4MIPs_2018-11-01/CMIP6/ScenarioMIP/UoM/UoM-AIM-ssp370-1-2-0/atmos/yr'
      $> target='./Out/trgas_ssp370_'

      and then I can call the routine like so:
      $> ./GHG_UKCA.py -s $source -o $target -p "scenarioMIP" -b 2015 -e 2035

      Alternatively, for the historical data:
      $> source='/gws/nopw/j04/cmip6_prep_vol1/cmip6_ancils/data/inputs4MIPs_2018-11-01/CMIP6/CMIP/UoM/UoM-CMIP-1-2-0/atmos/yr'
      $> target='./Out/trgas_historical_'

      and then I can call the routine like so:
      $> ./GHG_UKCA.py -s $source -o $target -p "historical" -b 1900 -e 2000

      *** NOTE! ***
      As of 25/11/18, the units attribute for the GHG variables is erroneously set to 1e-12 for all GHGs in the current source files. This has been raised with Paul Durack at LLNL. As a temporary fix, this script will write the correct units in the header of the output file without reading them from the source file. 

     T. Kuhlbrodt 21/11/18

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


# agreed format for columns that are not used by UKCA
ZERO = 0.0
DUM = -1.0e6

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

# temporary fix: write the units into the header of the output file
UNITSHEADER = "\n   UNITS:      ppt   ppm   ppb   ppb   ppt   ppt   ppt" \
    "ppt   ppt   ppt   ppt   ppt   ppt   ppt   ppt   ppt   ppt   ppt" \
    "   ppt   ppt   ppt   ppt   ppt   ppt"


def ukcaancil(start, end, source_files, output_file, project, data_version):
    # fixed parts of the source filename
    # this might change with the version of source files used
    PART1 = "mole_fraction_of_"
    PART2 = "_in_air"
    if project == "historical":
        if data_version is None:
            data_version = "v20160830"
    elif project == "scenarioMIP":
        if data_version is None:
            data_version = "v20180611"
    else:
        print("MIP project needs to be either 'historical' or 'scenarioMIP'.")
        exit()
    PART3 = os.path.join("gr1-GMNHSH", data_version)

    # calendar parameters
    if project == "historical":
        YRLEN = 365.24  # gregorian calendar
        YRBEG = 0
        YRSHIFT = 0
    elif project == "scenarioMIP":
        YRLEN = 365.0
        YRBEG = 1850
        YRSHIFT = 2015

    # Check that start is valid, else netcdf slice returns nothing (silently).
    # Output file will have corrected start year in filename, so the user is
    # informed of the change in two ways (filename and python warning).
    if start < YRSHIFT:
        warn(("'start' argument {} is earlier than origin year. Resetting "
              "to {}").format(start, YRSHIFT))
        start = YRSHIFT

    # load files
    molef = []
    units = []
    times = []
    for gas in GASES:
        factor = 1.0
        varname = PART1 + STDIC[gas][0] + PART2
        if project == "historical":  # turn underscores into dashes
            varname_mod = PART1 + STDIC[gas][0] + PART2
            varname_mod = varname_mod.replace("_", "-")
            fname = glob.glob(
                os.path.join(source_files[0], varname_mod, PART3, "*.nc"))
        elif project == "scenarioMIP":
            fname = glob.glob(
                os.path.join(source_files[0], varname, PART3, "*.nc"))
        print(fname[0])
        if len(fname) != 1:
            raise Exception("Either too many or no input nc files " +
                            gas["name"])
        #   add some exceptions
        if gas[0:4] == "DUMM":
            factor = DUM
        if gas == "HALON1202":
            factor = ZERO

        fh = Dataset(fname[0], mode="r")
        if end is None:
            molef.append(fh.variables[varname][start - YRSHIFT, 0] * factor)
            times.append(fh.variables["time"][start - YRSHIFT])
        else:
            molef.append(
                fh.variables[varname][start - YRSHIFT:end - YRSHIFT + 1, 0] *
                factor)
            times.append(fh.variables["time"][start - YRSHIFT:end - YRSHIFT +
                                              1])
        units.append(fh.variables[varname].units)
        fh.close()

    size = np.shape(molef)

    taxis = times[0]

    #   now write data to ASCII file
    if end is None:
        fname1 = output_file + str(start) + ".dat"
    else:
        fname1 = output_file + str(start) + "_" + str(end) + ".dat"
    with open(fname1, "w") as outp:
        #   the header
        outp.write("! " + project + " GHG data for CMIP6 \n! source: " +
                   source_files[0] + " \n")
        outp.write(" &THISFILE_SPECIFICATIONS \n")
        outp.write(" THISFILE_DATACOLUMNS    =          " + str(size[0]) +
                   ", \n")
        outp.write(" THISFILE_FIRSTYEAR      =        " + str(start) + ", \n")
        if end is None:
            outp.write(" THISFILE_LASTYEAR       =        " + str(start) +
                       ", \n")
        else:
            outp.write(" THISFILE_LASTYEAR       =        " + str(end) +
                       ", \n")
        outp.write(" THISFILE_ANNUALSTEPS    =           1, \n")
        outp.write(" THISFILE_FIRSTDATAROW   =          15, \n")
        outp.write(" THISFILE_UNITS          =  SEE ROW 13, \n")
        outp.write(" THISFILE_DATTYPE        =  RCPDAT \n/ \n")
        #   the data
        outp.write("   COLUMN:   ")
        for i in np.arange(size[0]):
            outp.write("   %i" % (i + 1), )
        #  temporary fix: write the units into the header of the output file
        outp.write(UNITSHEADER)
        #        outp.write("\n   UNITS:   ")
        #        for counter, value in enumerate(units):
        #            outp.write("   %s" %OMAG[value][0],)
        outp.write("\n   YEARS   ")
        for i in np.arange(size[0]):
            outp.write("   %s" % GASES[i], )
        if end is None:
            outp.write("\n   %i" % int(taxis / YRLEN + YRBEG), )
            for i in np.arange(size[0]):
                outp.write("    %7.4e" % molef[i])
        else:
            for j in np.arange(size[1]):
                outp.write("\n   %i" % int(taxis[j] / YRLEN + YRBEG), )
                for i in np.arange(size[0]):
                    outp.write("    %7.4e" % molef[i][j])
        outp.write("\n")


def _process(start, end, source_files, output_file, project, data_version):
    # Distinguish between timeslice/timeseries based on whether an end
    # argument is provided:
    ukcaancil(start, end, source_files, output_file, project, data_version)


if __name__ == "__main__":
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
