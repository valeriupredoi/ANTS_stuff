"""Tests for nitrogendeposition CLI.

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

from nitrogendeposition import (
    cmip6_ndep_jasmin_vn21,
    cmip6_ndep_jasmin_vn23_ssp585,
    cmip6_ndep_jasmin_vn25_past_and_future,
    cmip6_ndep_jasmin_vn26_past_and_future,
)

from nitrogendeposition.cmip6_ndep_jasmin_vn21 import main as vn21_main

# use these imports to test the executables (when/if they become executables)
# from nitrogendeposition.cmip6_ndep_jasmin_vn23_ssp585 import main as vn23_main
# from nitrogendeposition.cmip6_ndep_jasmin_vn25_past_and_future import main as vn25_main
# from nitrogendeposition.cmip6_ndep_jasmin_vn26_past_and_future import main as vn26_main
# from nitrogendeposition.cmip6_ndep import main as ndep_main


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


@patch('nitrogendeposition.cmip6_ndep_jasmin_vn21.main',
       new=wrapper(cmip6_ndep_jasmin_vn21))
def test_run_cmip6_ndep_jasmin_vn21_command():
    """Test run command."""
    executable = str(Path(cmip6_ndep_jasmin_vn21.__file__))
    print(f"Path to executable: {executable}")

    # test run as executable (executable functionality not yet set up!)
    # this can be set up in the main package via a list og console entries
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


def test_run_cmip6_ndep_jasmin_vn23_ssp585_command():
    """Test run command."""
    executable = str(Path(cmip6_ndep_jasmin_vn23_ssp585.__file__))
    print(f"Path to executable: {executable}")
    cmd = [executable, '--help']
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    stdoutdata, stderrdata = process.communicate()
    assert process.returncode == 0
    cmd = [executable, '--cow']
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    stdoutdata, stderrdata = process.communicate()
    assert process.returncode == 2


def test_run_cmip6_ndep_jasmin_vn25_past_and_future_command():
    """Test run command."""
    executable = str(Path(cmip6_ndep_jasmin_vn25_past_and_future.__file__))
    print(f"Path to executable: {executable}")
    cmd = [executable, '--help']
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    stdoutdata, stderrdata = process.communicate()
    assert process.returncode == 0
    cmd = [executable, '--cow']
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    stdoutdata, stderrdata = process.communicate()
    assert process.returncode == 2


def test_run_cmip6_ndep_jasmin_vn26_past_and_future_command():
    """Test run command."""
    executable = str(Path(cmip6_ndep_jasmin_vn26_past_and_future.__file__))
    print(f"Path to executable: {executable}")
    cmd = [executable, '--help']
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    stdoutdata, stderrdata = process.communicate()
    assert process.returncode == 0
    cmd = [executable, '--cow']
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    stdoutdata, stderrdata = process.communicate()
    assert process.returncode == 2
