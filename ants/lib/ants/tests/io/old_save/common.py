# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import subprocess
import unittest.mock as mock

from ants.deprecations import save_deprecation_message


def run_command(command, exception=None):
    """
    Run the command in a new process using :class:`subprocess.Popen`.

    Parameters
    ----------
    command : list of strings
        The command to run.
    exception : Exception
        The exception used if a non-zero return code is returned; if
        None, RuntimeError is used.

    Returns
    -------
    : str
        The standard output from the command.
    """
    with subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    ) as process:
        (stdoutdata, stderrdata) = process.communicate()
        return_code = process.returncode
        if return_code != 0:
            if exception is None:
                exception = RuntimeError
            raise exception(f"Return code: {return_code}: {stderrdata}")
    return stdoutdata


def ants_save_call():
    return mock.call(
        save_deprecation_message("ants.save"),
        FutureWarning,
        stacklevel=2,
    )


def ancil_save_call():
    return mock.call(
        save_deprecation_message("ants.fileformats.ancil.save"),
        FutureWarning,
        stacklevel=2,
    )


def netcdf_save_call():
    return mock.call(
        save_deprecation_message("ants.fileformats.netcdf.cf.save"),
        FutureWarning,
        stacklevel=2,
    )


def ukca_netcdf_save_call():
    return mock.call(
        save_deprecation_message("ants.fileformats.netcdf.ukca.save"),
        FutureWarning,
        stacklevel=2,
    )


def saver_call():
    return mock.call(
        save_deprecation_message("--saver"),
        FutureWarning,
        stacklevel=2,
    )
