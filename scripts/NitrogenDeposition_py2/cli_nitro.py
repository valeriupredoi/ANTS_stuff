"""
Module to provide a consistent command line interface.

Providing a consistent command line interface for many contrib applications
makes them easier to use, and easier to integrate into rose suites.

This should be imported into another module for use - see example.py for a
template.

CMIP6 specific instructions
===========================

Eventually, this will be catered for in the rose suite.  For now, though, two
things need to be done manually.  Firstly, check out the cli tool:

  https://code.metoffice.gov.uk/trac/ancil/browser/contrib/trunk/cmip6_utils/cli

Secondly, include the directory it was checked out into in your PYTHONPATH,
e.g.:
  export PYTHONPATH=/group_workspaces/jasmin2/tids/CMIP6_ANCIL/users/till/cli

"""

import argparse


def get_arg_parser(docs=None):
    """
    Provides an ArgumentParser with some pre-configured command line arguments.

    The parser can be extended with additional arguments.  After adding
    additional arguments (if any), simply call the parse_args method.  For
    example:

    arg_parser = get_arg_parser()
    args = arg_parser.parse_args()

    The provided arguments are:

    -s/--source  Source files for processing
    -o/--output  Output file to write the result to
    -n/--name    Species name for the preprocessing
    -b/--begin   Start year for the preprocesssing
    -e/--end     Final year for the preprocessing (optional)

    It is assumed that the processing will produce a timeslice if a "begin"
    argument is given without an "end" argument, and a timeseries if both
    "begin" and "end" are provided.  The "end" argument defaults to None so
    the processing code needs to be able to recognise that this means a
    timeslice should be generated.

    """

    # Creates ArgumentParser with output of the -h/--help argument set to the
    # docs argument:
    parser = argparse.ArgumentParser(
        conflict_handler="resolve",
        description=docs,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # The following two arguments will display in the group of "optional
    # arguments".  This is where arguments are added by default.

    # Add a single integer argument for the end year.  Also shows how to
    # define a default value for when the option isn't provided on the command
    # line.
    parser.add_argument(
        "-e",
        "--end",
        type=int,
        default=None,
        help=(
            "Optional end year for the processing.  If this is used, output "
            "is a timeseries from the 'begin' year to this year, inclusive "
            "of both the 'begin' year and the 'end' year.  Otherwise, the "
            "output will be a timeslice in the 'begin' year only."
        ),
    )

    parser.add_argument(
        "-n",
        "--name",
        type=str,
        help="Name of the species to be preprocessed.",
    )

    # We want to display some arguments as required, so we need to create a
    # separate group.  This group is used in the same way as the default group:
    required = parser.add_argument_group("required arguments")

    # Add a single integer argument for the start year.  Can be defined on
    # command line as either -b or --begin.  The --begin form defines the
    # variable in which the start year will be stored.
    #
    # More info on add_argument() can be found here:
    # https://docs.python.org/2/library/argparse.html#the-add-argument-method
    required.add_argument(
        "-b",
        "--begin",
        required=True,
        type=int,
        help="Start year for the processing.",
    )

    # Add a single string argument for the output file.
    required.add_argument(
        "-o",
        "--output",
        required=True,
        type=str,
        help="Filename to write the result to.",
    )

    # Add one or more string arguments for the input files.
    required.add_argument(
        "-s",
        "--sources",
        # SKL   required=True,
        required=False,
        nargs="+",
        type=str,
        help="Filename(s) to read the source data from.",
        dest="sources",
    )

    return parser
