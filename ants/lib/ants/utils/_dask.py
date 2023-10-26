# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import logging

import dask.array as da
import iris
import numpy as np

_LOGGER = logging.getLogger(__name__)


def copy(dask_array):
    """
    Perform true copy of a dask array.

    A copy of a dask array is returned where possible.  This is equivalent to
    if we were to call `dask.array.from_array` on the original array-like
    object it contains.  Since dask arrays are immutable, dask overrides dask
    array copy behaviour.  Currently limited to the following cases, otherwise
    return itself:

    - Dask array defined by a single function call within it's graph.

    Warning
    -------
    This is for ants core library usage ONLY!

    See Also
    --------
    https://code.metoffice.gov.uk/trac/ancil/ticket/1098 : for use case that
        requires this.

    """
    # TODO: Should no longer be needed after migration to iris 3 - see
    # e.g. #1673.
    if not is_lazy_data(dask_array):
        raise ValueError("array-like object appears not to be lazy.")
    res = dask_array
    array_obj = dask_array.dask.items()[-1][-1]
    _LOGGER.info(
        f"Dask array object: {array_obj}; "
        f"Dask dependencies: {dask_array.dask.dependencies}; "
        f"Dask chunksize: {dask_array.chunksize}"
    )
    if (
        len(dask_array.dask.dependencies) == 1
        and list(dask_array.dask.dependencies.values()) == [set()]
        and not isinstance(array_obj, tuple)
    ):
        res = da.from_array(array_obj, chunks=dask_array.chunksize, asarray=False)
    return res


def is_lazy_data(data):
    """
    Return whether the argument is a 'lazy' data array.

    At present, this means simply a Dask array.
    For now, this simply utilises the private iris function for this.

    Parameters
    ----------
    array : array-like
        An array-like object to query its lazy status.

    Returns
    -------
    : bool
        Returning True for lazy and False not.

    Warning
    -------
    This is for ants core library usage ONLY!

    """
    return iris._lazy_data.is_lazy_data(data)


def as_lazy_data(data, chunks=None, asarray=False):
    """
    Convert the input array `data` to a dask array.

    This is a placehlder for requesting lazy versions of an array-like objects.
    For now, this simply utilises the private iris function for this.

    Parameters
    ----------
    data : array-like
        An array. This will be converted to a dask array.
    chunks : int or (int, int, ...) or ((int, ...), (int, ...), ...), optional
        Describes how the created dask array should be split up. Defaults to a
        value first defined in biggus (being `8 * 1024 * 1024 * 2`).
        For more information see
        http://dask.pydata.org/en/latest/array-creation.html#chunks.
    asarray : bool
        If True, then chunks will be converted to instances of `ndarray`.
        Set to False (default) to pass passed chunks through unchanged.

    Returns
    -------
    : :class:`dask.array.array`
        The input array converted to a dask array.

    Warning
    -------
    This is for ants core library usage ONLY!

    """
    return iris._lazy_data.as_lazy_data(data, chunks=chunks, asarray=asarray)


def _is_masked_array(dask_array):
    """
    Return whether the chunktype of the input dask array is a np.ma.MaskedArray.

    This checks whether the _meta attribute of the given dask array is an
    instance of np.ma.MaskedArray.

    Parameters
    ----------
    dask_array : :class:`~dask.array.Array`
        A dask array to find out its chunktype.

    Returns
    -------
    : bool
        Returning True for a chunktype of np.ma.MaskedArray and False for anything
        else.

    Warning
    -------
    This is for ants core library usage ONLY! Behaviour of the _meta attribute is
    likely to change when the version of dask used is updated. See
    https://blog.dask.org/2019/06/22/dask-2.0 and
    https://docs.dask.org/en/latest/changelog.html#id195.

    """
    return isinstance(dask_array._meta, np.ma.MaskedArray)
