import os

import iris
import numpy as np
import pytest
from iris.cube import Cube, CubeList
from netCDF4 import Dataset

from greenhousegases.GHG_reform import interpol, main, write_rose_conf


def _create_sample_cube(varname, ancient_time=False):
    """Create sample cube."""
    data = np.ones((5, 10, 10), dtype=np.float32)
    cube = Cube(data, var_name=varname, units="J")
    if ancient_time:
        cube.add_dim_coord(
            iris.coords.DimCoord(
                [1000, 2000, 3000, 4000, 5000],
                standard_name="time",
                units="years since 01-01-01",
            ),
            0,
        )
    else:
        cube.add_dim_coord(
            iris.coords.DimCoord(
                [0, 1000, 2000, 3000, 4000],
                standard_name="time",
                units="years since 2000-01-01",
            ),
            0,
        )
    cube.add_dim_coord(
        iris.coords.DimCoord(
            [i + 0.5 for i in range(10)],
            standard_name="longitude",
            bounds=[[i, i + 1.0] for i in range(10)],
            units="degrees_east",
        ),
        1,
    )
    cube.add_dim_coord(
        iris.coords.DimCoord(
            [i + 0.5 for i in range(10)],
            standard_name="latitude",
            bounds=[[i, i + 1.0] for i in range(10)],
            units="degrees_north",
        ),
        2,
    )
    cube.attributes["STASH"] = "m01s00i030"

    return cube


# FIXME needs a test for all the globs


def test_interpolate():
    """Test interpolate() func."""
    array_in = np.ones((4))
    array_out = interpol(array_in)
    expected = np.ones((4))
    np.testing.assert_array_equal(expected, array_out)


def test_main(tmp_path):
    """Test ghg_ssp() func with project=historical."""
    # in data
    in_cube = _create_sample_cube("gasvar")
    in_cube = in_cube[0]
    gas_in_cube = _create_sample_cube("varname")
    save_path_root = tmp_path / "mole_fraction_of__in_air"
    save_path_gas = os.path.join(save_path_root, "HFC125")
    os.makedirs(save_path_gas)
    print("Test files location:", save_path_gas)
    cube_saved_path = os.path.join(save_path_gas, "indata.nc")
    print("Cube saved to:", cube_saved_path)
    iris.save(in_cube, cube_saved_path)

    start = 1
    end = 4

    gas_in = {
        "name": "HFC125",
        "converfac": 2,
        "units": 2,
        "varname": "gasvar",
        "sourcedir": save_path_root,
    }
    gas_out, time_out = main(gas_in)
    expected_gas_out = [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0]
    expected_time_out = [1849, 1850, 1851, 1852, 1853, 1854, 1855, 1856, 1857, 1858]
    np.testing.assert_array_equal(gas_out, expected_gas_out)
    np.testing.assert_array_equal(time_out, expected_time_out)


def test_write_rose_conf(tmp_path):
    """Test write_rose_conf() func."""
    # GASES = [CO2, CH4, N2O, CFC12, HFC]
    output_file = str(tmp_path)
    gas_mmr = {
        ("CO2", "mmr"): [3, 4, 5, 6, 7, 8, 9, 10],
        ("CO2", "year"): [1999, 2000, 2001, 2002, 2003, 2004],
        ("CH4", "mmr"): [3, 4, 5, 6, 7, 8, 9, 10],
        ("CH4", "year"): [1999, 2000, 2001, 2002, 2003, 2004],
        ("N20", "mmr"): [3, 4, 5, 6, 7, 8, 9, 10],
        ("N20", "year"): [1999, 2000, 2001, 2002, 2003, 2004],
        ("CFC12", "mmr"): [3, 4, 5, 6, 7, 8, 9, 10],
        ("CFC12", "year"): [1999, 2000, 2001, 2002, 2003, 2004],
        ("HFC134A", "mmr"): [3, 4, 5, 6, 7, 8, 9, 10],
        ("HFC134A", "year"): [1999, 2000, 2001, 2002, 2003, 2004],
        "sourcedir": output_file,
    }
    start = 1900
    end = 2200
    project = "historical"
    conf = write_rose_conf(gas_mmr)
    tested_file = str(tmp_path / "global_ghg_cmip6.conf")
    print("Output written to:", output_file)

    assert os.path.isfile(tested_file)
    with open(tested_file, "r") as file:
        all_lines = file.readlines()
        header = "[namelist:clmchfcg]\n"
        rando_string = "clim_fcg_nyears_n2o=6\n"
        assert header in all_lines
        assert rando_string in all_lines
