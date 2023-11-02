"""Tests for NitrogenDeposition CLI.

Includes a context manager to temporarily modify sys.argv
requirement: Python 3.7
"""
import subprocess
import contextlib
import copy
import functools
import sys

from pathlib import Path
from unittest.mock import patch
from io import StringIO

import pytest

from NitrogenDeposition import (
    cmip6_ndep_jasmin_vn21,
    cmip6_ndep_jasmin_vn23_ssp585,
    cmip6_ndep_jasmin_vn25_past_and_future,
    cmip6_ndep_jasmin_vn26_past_and_future,
    cmip6_ndep
)

from NitrogenDeposition.cmip6_ndep_jasmin_vn21 import main as vn21_main
from NitrogenDeposition import cmip6_ndep_jasmin_vn23_ssp585
from NitrogenDeposition import cmip6_ndep_jasmin_vn25_past_and_future
from NitrogenDeposition import cmip6_ndep_jasmin_vn26_past_and_future
from NitrogenDeposition import cmip6_ndep


def wrapper(f):
    @functools.wraps(f)
    def empty(*args, **kwargs):  # noqa
        if kwargs:
            raise ValueError(f'Parameters not supported: {kwargs}')
        return True

    return empty


@contextlib.contextmanager
def arguments(*args):
    backup = sys.argv
    sys.argv = list(args)
    print(backup)
    yield
    sys.argv = backup


def test_setargs():
    original = copy.deepcopy(sys.argv)
    with arguments('testing', 'working', 'with', 'sys.argv'):
        assert sys.argv == ['testing', 'working', 'with', 'sys.argv']
    assert sys.argv == original


@contextlib.contextmanager
def capture_sys_output():
    capture_out, capture_err = StringIO(), StringIO()
    current_out, current_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = capture_out, capture_err
        yield capture_out, capture_err
    finally:
        sys.stdout, sys.stderr = current_out, current_err


@patch('NitrogenDeposition.cmip6_ndep_jasmin_vn21.main', new=wrapper(cmip6_ndep_jasmin_vn21))
def test_run_cmip6_ndep_jasmin_vn21_command():
    """Test run command."""
    executable = str(Path(cmip6_ndep_jasmin_vn21.__file__))
    print(f"Path to executable: {executable}")

    # test run as executable (executable functionality not yet set up!)
    with arguments(executable, '--help'):
        cmip6_ndep_jasmin_vn21.main()
    with arguments(executable, '--cow'):
        cmip6_ndep_jasmin_vn21.main()

    # test run from command line
    cmd = [executable, '--help']
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    stdoutdata, stderrdata = process.communicate()
    assert process.returncode == 0
    cmd = [executable, '--cow']
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    stdoutdata, stderrdata = process.communicate()
    assert process.returncode == 2
