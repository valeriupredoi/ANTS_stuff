#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
"""
Mule cumf python utility
************************
This script is intended to provide access to basic mule cumf capability by
calling mule under the hood.
Currently this script exists to avoid the issues around having to build and
deploy um_utils and its dependencies under a given suite as utilising the
rose-ana MuleCumf class requires usage against the rose python environment
(that is python2.6).

"""
import argparse
import io
import os

import mule
from um_utils.cumf import COMPARISON_SETTINGS, UMFileComparison, full_report


def gen_cumf_report(source_filename, reference_filename):
    # Hardwire ignoring software version header entries as these are
    # valid differences
    flh_model_version_index = 12
    lookup_lbsrce_index = 38
    COMPARISON_SETTINGS["ignore_templates"] = {
        "fixed_length_header": [flh_model_version_index],
        "lookup": [lookup_lbsrce_index],
    }

    source_umfile = mule.AncilFile.from_file(source_filename)
    reference_umfile = mule.AncilFile.from_file(reference_filename)
    return UMFileComparison(source_umfile, reference_umfile)


def main(source_filename, reference_filename):

    assert os.path.isfile(source_filename), "No file to compare"
    assert os.path.isfile(reference_filename), "No file to compare against"

    comp = gen_cumf_report(source_filename, reference_filename)
    if not comp.match:
        output = io.StringIO()
        full_report(comp, stdout=output)
        contents = output.getvalue()
        output.close()
        msg = "FAIL: Ancillaries do not match: \n{}"
        raise RuntimeError(msg.format(contents))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    msg = "Source ancillary file."
    parser.add_argument("source_file", type=str, help=msg)

    msg = "Reference ancillary file"
    parser.add_argument("reference_file", type=str, help=msg)

    args = parser.parse_args()
    print(args)
    main(args.source_file, args.reference_file)
