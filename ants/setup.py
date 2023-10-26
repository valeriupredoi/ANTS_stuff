# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import importlib
import os
import warnings
from glob import glob

from setuptools import setup


def runnable_scripts(path):
    files = glob(os.path.join(path, "*"))
    scripts = [fnme for fnme in files if os.path.isfile(fnme)]
    return scripts


if __name__ == "__main__":
    setup(
        scripts=runnable_scripts("bin"),
    )

    for lib in ["gdal", "f90nml", "stratify", "ESMF"]:
        package = importlib.util.find_spec(lib)
        if not package:
            msg = 'Optional dependency "{}" not importable, see INSTALL.md'
            warnings.warn(msg.format(lib))
