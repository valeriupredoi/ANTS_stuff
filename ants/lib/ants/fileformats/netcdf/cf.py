# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import logging

import ants.utils
import iris.fileformats.netcdf as inetcdf
import numpy as np
from ants.config import CONFIG
from ants.deprecations import (
    issue_save_deprecation,
    save_deprecation_message_for_docstring,
)

_LOGGER = logging.getLogger(__name__)


def _coerce_netcdf_classic_dtypes(cubes):
    """
    Coerce cube data dtype to a netCDF classic compatible and CF compliant type.

    NetCDF classic types include:
        ['S1', 'i1', 'u1', 'i2', 'u2', 'i4', 'u4', 'i8', 'u8', 'f4', 'f8']
    We also ensure CF compliance by coercing unsigned types to signed types.

    Raises
    ------
    OverflowError
        If the array cannot be safely cast to the corresponding netCDF-classic
        type.

    .. warning::

        Coercing the data of the cube will load it into memory.

    """

    def _recast(array, target_type):
        if not np.can_cast(array.dtype, target_type):
            # Not generally cast safe so check whether we can coerce based on
            # the values actually present.
            maxval = array.max()
            minval = array.min()
            if (
                np.can_cast(maxval, target_type) is False
                or np.can_cast(minval, target_type) is False
            ):
                msg = (
                    "Cannot safely re-cast {} array to {} for writing to a "
                    "netCDF classic file"
                )
                raise OverflowError(msg.format(array.dtype, target_type))
        return array.astype(target_type)

    dconv = {
        "bool": "i1",  # byte
        "uint8": "i2",  # short
        "uint16": "i4",  # int
        "uint32": "i8",
        "uint64": "i8",
        "float16": "f4",
    }

    cubes = ants.utils.cube.as_cubelist(cubes)
    for cube in cubes:
        if cube.dtype.name in dconv:
            if (
                "valid_range" not in cube.attributes
                and "valid_min" not in cube.attributes
            ):
                if cube.dtype.kind == "b":
                    cube.attributes["valid_range"] = [0, 1]
                elif cube.dtype.kind == "u":
                    cube.attributes["valid_range"] = [0, np.iinfo(cube.dtype).max]

            target_type = dconv[cube.dtype.name]
            cube.data = _recast(cube.lazy_data(), target_type)

        for coord in cube.coords():
            if coord.dtype.name in dconv:
                target_type = dconv[coord.dtype.name]
                coord.points = _recast(coord.lazy_points(), target_type)
            if coord.has_bounds() and coord.bounds_dtype.name in dconv:
                target_type = dconv[coord.bounds_dtype.name]
                coord.bounds = _recast(coord.lazy_bounds(), target_type)


def _rechunk(cube):
    """
    Rechunks cube data in place such that innermost dimension is a single chunk.

    All other dimensions are automatically chunked.

    Parameters
    ----------
    cube : :class:`~iris.cube.Cube`

    Returns
    -------
    : None
    Operates on cube in place.
    """
    innermost_dimension = len(cube.shape) - 1
    rechunking = {i: "auto" for i in range(len(cube.shape))}
    rechunking[innermost_dimension] = -1
    cube.data = cube.core_data().rechunk(rechunking)


def _iris_netcdf4_classic_workaround(cubes):
    cubes = ants.utils.cube.as_cubelist(cubes)
    for cube in cubes:
        # If cubes are being saved in NETCDF4_CLASSIC format, have a dtype of
        # int64 and have lazy data, then the data needs to be realised or
        # converted to int32 to avoid a NetCDF4 RuntimeError. The approach
        # taken here is to realise the data, as this issue should only arise
        # for landseamasks and the conversion is already handled in iris for
        # non-lazy data. This issue is fixed in iris 3.2. At that point, this
        # workaround can be removed.  See
        # https://code.metoffice.gov.uk/trac/ancil/ticket/1527
        if cube.dtype == np.int64 and cube.has_lazy_data():
            cube.data


def _iris_dask_chunking_workaround(cubes):
    # If the innermost dimension is chunked, the netCDF save performance is
    # unacceptably slow.  We work around this by rechunking the innermost
    # dimension such that the entire dimension is a single chunk.  See
    # https://github.com/SciTools/iris/issues/4448 or #1555 for more details.
    # This behaviour can be disabled via configuration.
    cubes = ants.utils.cube.as_cubelist(cubes)
    for cube in cubes:
        # Default for CONFIG["ants_tuning"]["disable_rechunking"] is None, to
        # which we cannot apply the `.lower` string method.  So convert that
        # default to a string so we can be tolerant of users using "True",
        # "true", etc to disable rechunking:
        disable_rechunking_config = CONFIG["ants_tuning"]["disable_rechunking"]
        if disable_rechunking_config is None:
            disable_rechunking_config = "false"
        enable_rechunking = (
            disable_rechunking_config.lower() != "true" and cube.has_lazy_data()
        )
        _LOGGER.info(
            f"Rechunking is: {enable_rechunking} for cube {cube.name()}.  "
            f'Config setting was {CONFIG["ants_tuning"]["disable_rechunking"]}'
            ", which evaluated to an enable_rechunking of "
            f'{disable_rechunking_config.lower() != "true"}'
            f" and lazy data was {cube.has_lazy_data()}."
        )

        if enable_rechunking:
            # Can assume that core_data is a dask array here, since
            # enable_rechunking is false if the data is realised.
            _LOGGER.info(
                f"Before rechunking, cube {cube.name()} has chunks of "
                f"{cube.core_data().chunks}"
            )
            _rechunk(cube)
            _LOGGER.info(
                f"After rechunking, cube {cube.name()} has chunks of "
                f"{cube.core_data().chunks}"
            )


def save(cubes, filename, *args, **kwargs):
    """
    NetCDF saver.

    Note that if the cubes being saved have lazy data, this saver has two
    potential side effects:

    1. Any cubes of datatype int64 and being saved as NETCDF4_CLASSIC are realised.
    2. The data for each cube is rechunked such that the innermost dimension is
       a single chunk.  This can be disabled via configuration.

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
    :func:`iris.fileformats.netcdf.save` : for the underlying iris save operation.

    """
    # Issue a deprecation warning for the 'ants.fileformats.netcdf.cf.save'
    # interface.
    issue_save_deprecation("ants.fileformats.netcdf.cf.save")

    if "netcdf_format" not in kwargs:
        kwargs["netcdf_format"] = "NETCDF4_CLASSIC"
    if "unlimited_dimensions" not in kwargs:
        kwargs["unlimited_dimensions"] = []
    if kwargs["netcdf_format"] == "NETCDF4_CLASSIC":
        _coerce_netcdf_classic_dtypes(cubes)
        _iris_netcdf4_classic_workaround(cubes)

    _iris_dask_chunking_workaround(cubes)

    inetcdf.save(cubes, filename, *args, **kwargs)


save.__doc__ = (
    save_deprecation_message_for_docstring("ants.fileformats.netcdf.cf.save")
    + save.__doc__
)
