import os

import iris
import numpy as np
import pytest
from iris.cube import Cube, CubeList

from netCDF4 import Dataset
from greenhousegases.GHG_radiation import (
    WATERSHED,
    interpolate,
    ghg_ssp,
)


def _create_sample_cube(varname):
    """Create sample cube."""
    data = np.ones((5, 5, 5), dtype=np.float32)
    cube = Cube(data, var_name=varname, units="J")
    cube.add_dim_coord(
        iris.coords.DimCoord(
            [0, 1, 2, 3, 4],
            standard_name="time",
            units="days since 2000-01-01",
        ),
        0,
    )
    cube.add_dim_coord(
        iris.coords.DimCoord(
            [i + 0.5 for i in range(5)],
            standard_name="longitude",
            bounds=[[i, i + 1.0] for i in range(5)],
            units="degrees_east",
        ),
        1,
    )
    cube.add_dim_coord(
        iris.coords.DimCoord(
            [i + 0.5 for i in range(5)],
            standard_name="latitude",
            bounds=[[i, i + 1.0] for i in range(5)],
            units="degrees_north",
        ),
        2,
    )
    cube.attributes["STASH"] = "m01s00i030"

    return cube


# FIXME needs a test for all the globs


def test_interpolate():
    """Test interpolate() func."""
    array_in = np.ones((2, 3, 4))
    array_out = interpolate(array_in)
    expected = np.ones((1, 3, 4))
    np.testing.assert_array_equal(expected, array_out)


def test_ghg_ssp_historical(tmp_path):
    """Test ghg_ssp() func with project=historical."""
    # in data
    in_cube = _create_sample_cube("mole_fraction_of__in_air")
    gas_in_cube = _create_sample_cube("varname")
    gas_in = {
        "name": "HFC125", "converfac": 2, "units": 2, "varname": ""
    }

    # run routine
    project = "historical"
    start = WATERSHED
    end = start + 2

    save_path = tmp_path / "mole-fraction-of--in-air"
    save_path = os.path.join(save_path, "gr1-GMNHSH/v20160830")
    os.makedirs(save_path)
    print("Test files location:", save_path)
    cube_saved_path = os.path.join(save_path, "indata.nc")
    print("Cube saved to:", cube_saved_path)
    iris.save(in_cube, cube_saved_path)
    source_files = str(tmp_path)

    with pytest.raises(SystemExit) as exc:
        gas_out, time_out = ghg_ssp(
            gas_in, start, end,
            source_files, project,
            data_version="v20160830"
        )
        srexc = "Starting year needs to be earlier than finishing year"
        assert srexc in str(exc)

    
