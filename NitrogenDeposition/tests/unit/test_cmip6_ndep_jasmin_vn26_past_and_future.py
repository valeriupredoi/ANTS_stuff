import os
import sys

import subprocess
import contextlib
import functools
import iris
import numpy as np
import pytest

from pathlib import Path
from unittest.mock import patch
from iris.cube import Cube, CubeList
from nitrogendeposition import cmip6_ndep_jasmin_vn26_past_and_future
from nitrogendeposition.cmip6_ndep_jasmin_vn26_past_and_future import (regrid,
    load_cube_clim_noclim,
    save_cube,
    main
)


@contextlib.contextmanager
def arguments(*args):
    backup = sys.argv
    sys.argv = list(args)
    yield
    sys.argv = backup


def wrapper(f):
    @functools.wraps(f)
    def empty(*args, **kwargs):
        if kwargs:
            raise ValueError(f'Parameters not supported: {kwargs}')
        return True

    return empty


def _create_sample_cube():
    """Create sample cube."""
    data = np.ones((2, 5, 5), dtype=np.float32)
    cube = Cube(data, var_name='co2', units='J')
    cube.add_dim_coord(
        iris.coords.DimCoord(
            [0, 1],
            standard_name='time',
            units='days since 2000-01-01',
        ),
        0,
    )
    cube.add_dim_coord(
        iris.coords.DimCoord(
            [i + .5 for i in range(5)],
            standard_name='longitude',
            bounds=[[i, i + 1.] for i in range(5)],
            units='degrees_east',
        ),
        1,
    )
    cube.add_dim_coord(
        iris.coords.DimCoord(
            [i + .5 for i in range(5)],
            standard_name='latitude',
            bounds=[[i, i + 1.] for i in range(5)],
            units='degrees_north',
        ),
        2,
    )
    cube.attributes["STASH"] = "m01s00i030"

    return cube


def test_regrid(tmp_path):
    """Test regrid."""
    in_cube = _create_sample_cube()
    grid_save_path = tmp_path / "gridFile.nc"
    iris.save(in_cube, str(grid_save_path))
    tested_out = regrid(in_cube, str(grid_save_path))
    np.testing.assert_array_equal(tested_out.data, in_cube.data)


def test_load_cube_clim_noclim(tmp_path):
    """Test load_cube_clim_noclim."""
    in_cube = _create_sample_cube()
    grid_save_path = tmp_path / "gridFile.nc"
    iris.save(in_cube, str(grid_save_path))
    filenames = [str(grid_save_path)]
    tested_out = load_cube_clim_noclim(filenames, clim=False)
    np.testing.assert_array_equal(tested_out[0].data, in_cube.data)
    tested_out = load_cube_clim_noclim(filenames, clim=True)
    assert not tested_out


@patch('nitrogendeposition.cmip6_ndep_jasmin_vn26_past_and_future.main',
       new=wrapper(main))
def test_main():
    """Test main."""
    tested_out = main()
    assert not tested_out
    with arguments('python', 'nitrogendeposition/cmip6_ndep_jasmin_vn26_past_and_future.py',
                   '--end', '2000'):
        assert not main()


def test_save_cube(tmp_path):
    """Test save_cube."""
    in_cube = _create_sample_cube()
    filename = str(tmp_path / "oldsaver.nc")
    save_cube(in_cube, filename, use_new_saver=False)
    assert os.path.isfile(filename)
    filename = str(tmp_path / "newsaver.anc")
    save_cube(in_cube, filename, use_new_saver=True)
    assert os.path.isfile(filename)


def test_main_no_resolution():
    """Test various components of main."""
    executable = str(Path(cmip6_ndep_jasmin_vn26_past_and_future.__file__))
    print(f"Path to executable: {executable}")

    cmd = [executable, "--end", "2000"]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdoutdata, stderrdata = process.communicate()
    err = "the following arguments are required: -r/--resolution"
    assert err in str(stderrdata)


def test_main_withNdep():
    """Test various components of main."""
    executable = str(Path(cmip6_ndep_jasmin_vn26_past_and_future.__file__))
    print(f"Path to executable: {executable}")

    cmd = [
        executable,
        "--resolution", "./TESTDIR_N96E_SSP585/qrparm.mask_n96e_eorca1_v2.2x_ESMF",
        "-n", "Ndep"
    ]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdoutdata, stderrdata = process.communicate()
    print(stderrdata)
    err = "unrecognized arguments: -n Ndep"
    assert err in str(stderrdata)
