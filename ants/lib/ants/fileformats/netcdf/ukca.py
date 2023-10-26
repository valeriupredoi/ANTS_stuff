# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import warnings

import ants
import dask
import numpy as np
from ants.deprecations import (
    issue_save_deprecation,
    save_deprecation_message_for_docstring,
)

from . import cf

LOCAL_ATTS = [
    "tracer_name",
    "vertical_scaling",
    "lowest_level",
    "highest_level",
    "hourly_scaling",
]


def _update_conventions(cubes):
    """
    Update old UKCA conventions used by the cubes.

    This is only needed where we are regridding older master files which use
    the previous conventions. The convention changes are:
        1. Replace emission_type with update_type
        2. Encode integers as netcdf ints (32-bit) instead of strings

    """
    for cube in cubes:
        # Replace emission type with update type attribute.
        if "emission_type" in cube.attributes:
            if "update_type" not in cube.attributes:
                cube.attributes["update_type"] = cube.attributes["emission_type"]
            del cube.attributes["emission_type"]

        # convert all UKCA numeric attributes to 32-bit integers.
        int_attrib = [
            "update_type",
            "update_freq_in_hours",
            "lowest_level",
            "highest_level",
        ]
        for attrib in int_attrib:
            if attrib in cube.attributes:
                cube_attrib = cube.attributes[attrib]
                cube.attributes[attrib] = np.int32(cube_attrib)


def _ukca_conventions(cubes):
    cubes = ants.utils.cube.as_cubelist(cubes)
    _update_conventions(cubes)
    for cube in cubes:
        if cube.dtype != np.int32:
            cube.data = cube.core_data().astype(np.float32, copy=False)
        # If the data hasn't been realised, check whether the array contained in the
        # dask.array is a np.ma.MaskedArray to avoid realising the data.
        if cube.has_lazy_data():
            if ants.utils._dask._is_masked_array(cube.core_data()):
                warnings.warn(
                    "Cube has masked points. Filling with zeros as per UKCA convention."
                )
                cube.data = dask.array.ma.filled(cube.core_data(), fill_value=0)
        # If the data is already realised, can directly check the data for a mask.
        else:
            if np.ma.isMaskedArray(cube.data):
                warnings.warn(
                    "Cube has masked points. Filling with zeros as per UKCA convention."
                )
                cube.data = cube.data.filled(fill_value=0)


def save(cubes, filename, *args, **kwargs):
    """
    UKCA NetCDF saver.

    The following represents applied UKCA specification or standard UKCA
    NetCDF setup:

    - Compression used: zlib with complevel=4
    - Ensure that specific attributes are local (see :obj:`LOCAL_ATTS`).
    - Bounds are present (guessed where not).
    - Data types are made either 32-bit integer or 32-bit float.
    - Masked data is filled with zeros (a warning is issued where this
      happens).
    - Old UKCA conventions present are updated:
        - Replace `emission_type` with `update_type`.
        - All UKCA numeric attributes are converted to 32-bit integers.  These
          include `update_type`, `update_freq_in_hours`, `lowest_level` and
          `highest_level`.

    Parameters
    ----------
    cubes : One or more :class:`~iris.cube.Cube`
        Input cubes to save.
    filename : str
        Output filename.
    netcdf_format : :obj:`str`, optional
        Underlying netCDF file format, one of 'NETCDF4', 'NETCDF4_CLASSIC',
        'NETCDF3_CLASSIC' or 'NETCDF3_64BIT'. Default is 'NETCDF4_CLASSIC'
        format.
    unlimited_dimensions : iterable of str and/or :class:`iris.coords.Coord`
        Coordinate dimensions of `cube` to save with the NetCDF dimension
        variable length 'UNLIMITED'.  By default there are no dimensions
        assigned with length 'UNLIMITED'.

    See Also
    --------
    :func:`ants.fileformats.netcdf.cf.save` : for the generic netCDF save function.

    """
    # Issue a deprecation warning for the 'ants.fileformats.netcdf.ukca.save'
    # interface.
    issue_save_deprecation("ants.fileformats.netcdf.ukca.save")

    if "zlib" not in kwargs:
        kwargs["zlib"] = True
    if "complevel" not in kwargs:
        kwargs["complevel"] = 4
    if "local_keys" not in kwargs:
        kwargs["local_keys"] = LOCAL_ATTS
    else:
        kwargs["local_keys"] = set(list(kwargs["local_keys"]) + LOCAL_ATTS)

    ants.utils.cube.guess_horizontal_bounds(cubes)
    _ukca_conventions(cubes)
    cf.save(cubes, filename, *args, **kwargs)


save.__doc__ = (
    save_deprecation_message_for_docstring("ants.fileformats.netcdf.ukca.save")
    + save.__doc__
)
