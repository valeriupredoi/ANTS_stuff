#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
"""
nccmp comparison utility
*************************
This script wraps the nccmp commandline utility to compare many ancillaries in
a single call.

"""
import argparse
import os
import re
import subprocess

import iris


# This may still be required as there may still be situations where switching
# on the buffer flag is not a good idea.
def use_buffer_flag(filename):
    """
    Determine whether we should use the buffer (-B) flag or not with nccmp.

    The logical behaviour is as follows:
    - Return False if there are any unlimited dimensions (it has been shown
      that use of the buffer keyword isn't necessary under at least some
      circumstances where unlimited dimensions are specified.
    - Any large outer dimensions (large being >= 500) for > 2D datasets will
      return True for buffer flag usage.
    - Above rules not hit, so return False as our default (don't use the
      buffer keyword).

    See Also
    --------
    https://gitlab.com/remikz/nccmp/-/issues/5 : for nccmp issue.

    """
    dim_large_threshold = 500

    icffh = iris.fileformats.cf.CFReader(filename)

    # Iterate over each of our data variables so see if ANY qualify for a
    # decision as to buffer keyword to nccmp usage is appropriate.
    for var in icffh.cf_group.data_variables.items():
        # Any unlimited dimensions? - don't use buffer flag.
        if re.search("unlimited dimensions: .+\n", str(var)):
            return False

    for var in icffh.cf_group.data_variables.items():
        # Any high dimensional outer dimensions - use buffer flag.
        shape = eval(re.search("current shape = (.*)", str(var)).groups()[0])
        if len(shape) > 2 and any([shp > dim_large_threshold for shp in shape[:-2]]):
            return True

    # No conditions met above, then just don't use the buffer
    return False


def main(source_filename, reference_filename, exclude_attr):

    assert os.path.isfile(source_filename), "No file to compare"
    assert os.path.isfile(reference_filename), "No file to compare against"

    command = (
        "nccmp --var-diff-count 5 --global --data --force --metadata "
        "--nans-are-equal --buffer-entire-var --statistics"
    )
    command = f"{command} {source_filename} {reference_filename}"
    if exclude_attr:
        # exclude user defined global attributes and all local
        # history attributes
        command = "{} --globalex={}".format(command, ",".join(exclude_attr))

    # We cannot use subprocess.call as nccmp doesn't always return an exit code
    # we might expect. Instead, we check whether it returns any output at all.
    with subprocess.Popen(command, stderr=subprocess.PIPE, shell=True) as proc:
        diffs = proc.stderr.read()
    if diffs:
        msg = "NetCDF comparison failure.\nFile1:{}\nFile2:{}\n{}\n{}"
        raise RuntimeError(
            msg.format(source_filename, reference_filename, command, diffs)
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    msg = "Source ancillary file."
    parser.add_argument("source_file", type=str, help=msg)

    msg = "Reference ancillary file"
    parser.add_argument("reference_file", type=str, help=msg)

    msg = "Attributes to ignore"
    parser.add_argument("--exclude-attr", type=str, nargs="+", help=msg)
    args = parser.parse_args()
    main(args.source_file, args.reference_file, args.exclude_attr)
