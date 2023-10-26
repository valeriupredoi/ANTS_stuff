# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import numpy as np
import scipy.sparse.linalg
from scipy import sparse

# Stop rounding problems (close to the pole), we want to limit eps to be some
# fraction of the value at which eps + 1 is rounded to eps.
# This is currently set to be a 50th the maximum (32-bit) value of eps.
# This is to a great extent to remain consistent with CAP behaviour.
MAXEPS = 1 / (50.0 * np.finfo("float32").eps)


def _matrix_interior(size, epsilon):
    # https://doi.org/10.1175/1520-0493(1988)116<2132:HOLPIT>2.0.CO;2
    # Eqn. 15-18
    # Inside points
    z = np.ones(size) * (1.0 - epsilon)
    a = np.ones(size) * (6.0 * (1.0 + epsilon))
    b = np.ones(size) * (15.0 * (1.0 - epsilon))
    c = np.ones(size) * (20.0 * (1.0 + epsilon))
    d = b
    e = a
    f = z
    return z, a, b, c, d, e, f


def define_matrix_periodic(size, epsilon):
    """
    Return matrix 'C' in the linear system `Cphi = g`

    Parameters
    ----------
    size : size
        Length of 1D array to which the filter is to be applied.
    epsilon: float
        Epsilon determines at what scale the transition occurs.

    Returns
    -------
    : ndarray
        Return 'C' in the linear system `Cphi = g` where 'C' is a sparse csc
        array of shape when expanded (size x size).

    Notes
    -----
    Check Appendix A "Matrix Inversion-Cyclic Case" from the following:
    https://journals.ametsoc.org/doi/pdf/10.1175/1520-0493%281991% \
    29119%3C0477%3AARORAI%3E2.0.CO%3B2


    """
    z, a, b, c, d, e, f = _matrix_interior(size, epsilon)
    # Three extra diagonals.
    diagonals = [
        d[3:4],
        e[3:5],
        f[3:6],
        z[3:],
        a[2:],
        b[1:],
        c,
        d[:-1],
        e[:-2],
        f[:-3],
        z[3:6],
        a[3:5],
        b[3:4],
    ]
    offsets = np.array(
        [
            1 - size,
            2 - size,
            3 - size,
            -3,
            -2,
            -1,
            0,
            1,
            2,
            3,
            size - 3,
            size - 2,
            size - 1,
        ]
    )
    return sparse.diags(diagonals, offsets, shape=(size, size), format="csc")


def define_matrix_non_periodic(size, epsilon):
    """
    Return 'C' in the linear system `Cphi = g`

    Parameters
    ----------
    size : size
        Length of 1D array to which the filter is to be applied.
    epsilon: float
        Epsilon determines at what scale the transition occurs.

    Returns
    -------
    : ndarray
        Return 'C' in the linear system `Cphi = g` where 'C' is a sparse csc
        array of shape when expanded (size x size).

    """
    # Variable names here are to be consistent with both paper and the CAP.
    z, a, b, c, d, e, f = _matrix_interior(size, epsilon)

    # Set boundary points.
    z[0:3] = z[:-4:-1] = 0.0
    a[0] = a[1] = 0.0
    a[2] = 1.0 + epsilon
    a[:-4:-1] = a[:3]
    b[0] = 0
    b[1] = 1.0 - epsilon
    b[2] = 4 * (1.0 - epsilon)
    b[:-4:-1] = b[:3]
    c[0] = 1
    c[1] = 2 * (1.0 + epsilon)
    c[2] = 6 * (1.0 + epsilon)
    c[:-4:-1] = c[:3]

    diagonals = [z[3:], a[2:], b[1:], c, d[:-1], e[:-2], f[:-3]]
    offsets = np.array([-3, -2, -1, 0, 1, 2, 3])
    return sparse.diags(diagonals, offsets, shape=(size, size), format="csc")


def filters_periodic(data, epsilon):
    """
    Return 'g' in the linear system `Cphi = g`

    6th-order periodic filter.

    Parameters
    ----------
    data : ndarray
        Data on which to apply this filter.
    epsilon: float
        Epsilon determines at what scale the transition occurs.

    Returns
    -------
    : ndarray
        Return 'g' in the linear system `Cphi = g` which is the same shape as
        the provided data array.

    """
    rhs = data.copy()
    datac = np.concatenate([data, data[0:3]], 0)
    datac = np.concatenate([data[-3:], datac], 0)
    rhs[:] = epsilon * (
        (datac[0:-6] + datac[6:])
        - 6.0 * (datac[1:-5] + datac[5:-1])
        + 15.0 * (datac[2:-4] + datac[4:-2])
        - 20.0 * datac[3:-3]
    )
    return rhs


def filters_non_periodic(data, epsilon):
    """
    Return 'g' in the linear system `Cphi = g`

    6th-order filter.

    Parameters
    ----------
    data : ndarray
        Data on which to apply this filter.
    epsilon: float
        Epsilon determines at what scale the transition occurs.

    Returns
    -------
    : ndarray
        Return 'g' in the linear system `Cphi = g` which is the same shape as
        the provided data array.

    """
    rhs = data.copy()
    rhs[3:-3] = epsilon * (
        (data[0:-6] + data[6:])
        - 6.0 * (data[1:-5] + data[5:-1])
        + 15.0 * (data[2:-4] + data[4:-2])
        - 20.0 * data[3:-3]
    )
    # Boundary filter conditions
    rhs[0] = rhs[-1] = 0

    rhs[1] = epsilon * (data[0] - 2 * data[1] + data[2])
    rhs[-2] = epsilon * (data[-3] - 2 * data[-2] + data[-1])

    rhs[2] = epsilon * (
        -1.0 * (data[0] + data[4]) + 4 * (data[1] + data[3]) - 6 * (data[2])
    )
    rhs[-3] = epsilon * (
        -1.0 * (data[-1] + data[-5]) + 4 * (data[-2] + data[-4]) - 6 * (data[-3])
    )
    return rhs


def epsilon_iso(shape, epsilon):
    """
    We calculate epsilon as a function of latitude and capture this as it's
    own function rather that through pure cython to aid testing.

    """
    # Variable names derived from the CAP.

    # Parameters used to derive epsilon to ensure isotropy
    # ----------------------------------------------------
    thetac = 2.0 * np.arctan(epsilon ** (-1.0 / 6.0))
    # xlenc: No. of grid-lengths at which power=0.5
    xlenc = 2.0 * np.pi / thetac
    # nsvew: resolution in north-south direction divided by that in east-west
    # direction
    nsvew = np.abs(2.0 * (shape[0]) / (shape[1]))

    # Derive epsilon
    inter = nsvew * np.abs(
        np.sin(((np.arange(1, shape[0] + 1)) - 0.5) * np.pi / (shape[0]))
    )
    # coslat: number of grid-boxes in the east-west direction that
    # span the same real distance as one grid-box in the
    # north-south direction.
    coslat = 1.0 / inter
    coslat[coslat < 1] = 1.0
    # wavelength: length-scale in grid-lengths at which we want 50%
    # filtering.
    wavelength = coslat * xlenc
    eps = 1.0 / (np.tan(np.pi / wavelength) ** 6)
    # Ensure isotropy applied only polewards of where gridboxes are
    # square.
    eps[inter == 0] = epsilon
    eps = np.minimum(eps, MAXEPS)

    return eps


def _fls_xlenc(filter_length_scale, delta_lambda, earth_radius):
    return 180.0 * filter_length_scale / (np.pi * earth_radius * delta_lambda)


def filter_length_scale_to_eps(filter_length_scale, delta_lambda, earth_radius):
    """
    Convert the filter length scale to epsilon in the north south direction.

    """
    # Variable names derived from the CAP.

    xlenc = _fls_xlenc(filter_length_scale, delta_lambda, earth_radius)
    edubl = 1.0 / ((np.tan(np.pi / xlenc)) ** 6)
    return edubl


def filter_length_scale_to_eps_iso(
    filter_length_scale, lats, delta_lambda, earth_radius
):
    """Calculate epsilon as a function of latitude."""
    # Variable names derived from the CAP.
    delta_phi = lats[1] - lats[0]

    xlenc = _fls_xlenc(filter_length_scale, delta_lambda, earth_radius)
    edubl = filter_length_scale_to_eps(filter_length_scale, delta_lambda, earth_radius)
    nsvew = delta_lambda / delta_phi

    coslat = 1.0 / (nsvew * abs(np.cos(np.deg2rad(lats))))
    wavelength = coslat * xlenc
    eps = 1.0 / ((np.tan(np.pi / wavelength)) ** 6)
    eps = np.minimum(eps, MAXEPS)
    eps[np.abs(lats) == 90] = edubl
    return eps


def filter_rows_varying_epsilon(data, eps, periodic):
    """
    Filter rows, utilising a varying epsilon value.

    Changes are made in place to the data.
    """
    for i in range(data.shape[0]):
        dat = data[i, :].reshape(-1)
        epsilon = eps[i]
        if periodic:
            filt_g = filters_periodic(dat, epsilon)
            arr_C = define_matrix_periodic(data.shape[1], epsilon)
        else:
            filt_g = filters_non_periodic(dat, epsilon)
            arr_C = define_matrix_non_periodic(data.shape[1], epsilon)
        perms_phi = scipy.sparse.linalg.spsolve(arr_C, filt_g)
        msg = "The filter (phi) shape {} does not match the data row shape {}."
        assert dat.shape == perms_phi.shape, msg.format(perms_phi.shape, dat.shape)
        dat[:] += perms_phi


def filter_rows_constant_epsilon(data, epsilon, periodic):
    """
    Filter rows, utilising a constant epsilon value.

    Changes are made in place to the data.

    Parameters
    ----------
    data : ndarray
        1D or 2D array to be filtered.
    epsilon: float
        Filter parameter epsilon to be used for each row.
    periodic : bool, optional
        Specify whether the array provided is to be handled as periodic or not.

    Returns
    -------
        In-place operation.

    """
    if periodic:
        arr_C = define_matrix_periodic(data.shape[1], epsilon)
        filt_g = filters_periodic(data.T, epsilon)
    else:
        arr_C = define_matrix_non_periodic(data.shape[1], epsilon)
        filt_g = filters_non_periodic(data.T, epsilon)
    perms_phi = scipy.sparse.linalg.spsolve(arr_C, filt_g)
    data[:] += perms_phi.T


def filter_columns(data, epsilon):
    """
    Filter columns, utilising a constant epsilon value.

    Changes are made in place to the data.

    Parameters
    ----------
    data : ndarray
        1D or 2D array to be filtered.
    epsilon: float
        Filter parameter epsilon to be used for each column.

    Returns
    -------
        In-place operation.

    """
    filt_g = filters_non_periodic(data, epsilon)
    arr_C = define_matrix_non_periodic(data.shape[0], epsilon)
    perms = scipy.sparse.linalg.spsolve(arr_C, filt_g)
    data[:] += perms


def raymond_filter_ndarray(data, eps_col, eps_row, periodic=False):
    """
    Raymond filter applied to a provided data array.

    Parameters
    ----------
    data : ndarray
        1D or 2D array to be filtered.
    eps_col : float
        Filter parameter epsilon to be used for each column.
    eps_row : float or iterale of floats
        Filter parameter(s) epsilon to be used for each row.
    periodic : bool, optional
        Specify whether the array provided is to be handled as periodic or not.

    Returns
    -------
        In-place operation.

    """
    if data.ndim == 1:
        data = data.reshape(1, -1)

    # Apply the filter to the data rows and columns.
    if data.shape[1] > 1 and np.asarray(eps_row).size < 2:
        filter_rows_constant_epsilon(data, eps_row, periodic)
    else:
        # Filter row by row as there is a different value of epsilon for each
        # row.
        filter_rows_varying_epsilon(data, eps_row, periodic)
    if data.shape[0] > 1:
        filter_columns(data, eps_col)
