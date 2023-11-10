"""Tests for NitrogenDeposition CLI.

Includes a context manager to temporarily modify sys.argv
requirement: Python 3.7
"""
import subprocess
# import contextlib
# import copy
# import functools
# import sys

from pathlib import Path
# from unittest.mock import patch
# from io import StringIO

import pytest

from greenhousegases import GHG_radiation
from greenhousegases import ghg_ukca_ncplot
from greenhousegases import GHG_UKCA


def test_run_GHG_radiation_command():
    """Test run command."""
    executable = str(Path(GHG_radiation.__file__))
    print(f"Path to executable: {executable}")
    cmd = [executable, '--help']
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    stdoutdata, stderrdata = process.communicate()
    assert process.returncode == 0
    cmd = [executable, '--cow']
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    stdoutdata, stderrdata = process.communicate()
    assert process.returncode == 2


@pytest.mark.skip(
    reason='Cannot run executable without access to data.'
)
def test_run_ghg_ukca_ncplot_command():
    """Test run command."""
    executable = str(Path(ghg_ukca_ncplot.__file__))
    print(f"Path to executable: {executable}")
    cmd = [executable, '--help']
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    stdoutdata, stderrdata = process.communicate()
    assert process.returncode == 0
    cmd = [executable, '--cow']
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    stdoutdata, stderrdata = process.communicate()
    assert process.returncode == 2


def test_run_GHG_UKCA_command():
    """Test run command."""
    executable = str(Path(GHG_UKCA.__file__))
    print(f"Path to executable: {executable}")
    cmd = [executable, '--help']
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    stdoutdata, stderrdata = process.communicate()
    assert process.returncode == 0
    cmd = [executable, '--cow']
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    stdoutdata, stderrdata = process.communicate()
    assert process.returncode == 2
