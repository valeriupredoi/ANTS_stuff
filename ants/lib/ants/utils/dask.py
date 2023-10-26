# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants
import dask.array as da


def deferred_data_update(data, newdata, slices):
    """
    Assign new data to the specified slices of data lazily.

    With the benefits of dask, override 'data' with 'newdata' at the
    specified slices without realising any data.

    Parameters
    ----------
    data : Array-like object
        2D array like object which includes numpy arrays or dask arrays.
        This represents the arrays which the other is transplanted onto.
    newdata : Array-like object
        2D array like object which includes numpy arrays or dask arrays.
        This represents the array which is transplanted into the other.
    slices : tuple(slice, slice)
        Slices object representing the 2D slicing of 'data' to transplant the
        'newdata'.

    Returns
    -------
    : :class:`dask.array`
        Lazy array representing the original array with the new data
        transplanted onto it.

    Note
    ----
    This function serves the purpose of replacing part of a dask array with
    some changes.

    """
    if data.ndim != 2:
        msg = "Expected 2D source data, got {} instead".format(data.ndim)
        raise ValueError(msg)
    if newdata.ndim != 2:
        msg = "Expected 2D target data, got {} instead".format(newdata.ndim)
        raise ValueError(msg)

    data = ants.utils._dask.as_lazy_data(data)
    idata = ants.utils._dask.as_lazy_data(newdata)

    ndata = ants.utils._dask.copy(data)[slice(slices[0].stop, None), slices[1]]
    if 0 not in ndata.shape:
        idata = da.concatenate([idata, ndata], 0)

    ndata = ants.utils._dask.copy(data)[
        slice(slices[0].start, None), slice(slices[1].stop, None)
    ]
    if 0 not in ndata.shape:
        idata = da.concatenate([idata, ndata], 1)

    ndata = ants.utils._dask.copy(data)[
        slice(slices[0].start, None), slice(None, slices[1].start)
    ]
    if 0 not in ndata.shape:
        idata = da.concatenate([ndata, idata], 1)

    ndata = ants.utils._dask.copy(data)[slice(None, slices[0].start), :]
    if 0 not in ndata.shape:
        idata = da.concatenate([ndata, idata], 0)
    return idata
