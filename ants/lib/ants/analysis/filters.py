# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants
import numpy as np

from . import _raymond


def raymond(source, epsilon=None, filter_length_scale=None, isotropic=False):
    """
    Raymond filter applied to a provided cube.

    This function is experimental and may be removed without notice.
    An in-place high-order low-pass implicit tangent 1D filter, applied to
    n-dimensions.

    This is achieved by solving a linear system (Ax=b), following
    common terminology for Raymond Filtering, `Cphi = g`, where 'g' corresponds
    to the filter (as returned by `filters`), 'C' corresponds to the matrix.
    'phi' corresponds to the set of permutations being solved by the filter
    and applied to the supplied cube.
    A latitude varying epsilon ensures an isotropic filtering is applied.  That
    is, since the grid spacing increases as we go away from the equator,
    adjustments are made to epsilon to ensure the filter is applied relative to
    this difference in scale.
    Periodic/cyclic fields are also supported (that is longitudinal
    wraparound).

    Where the filter is specifically requested to be isotropic, adjustments to
    epsilon/filter_length_scale are made away from the equator to ensure
    isotropic filter behaviour.  This is done in a way consistent with the CAP.
    This is currently limited to the source fields with the following
    characteristics::

      - Global
      - Regularly spaced.
      - Standard lat-lon coordinate system compatible with the UMSPHERE
        (see :func:`ants.coord_systems.as_ants_crs`).

    The tolerance for isotropy checks on the source fields can be configured
    (see :class:`ants.config.GlobalConfiguration`).

    Parameters
    ----------
    source : :class:`iris.cube.Cube`
        Field to have the raymond filter applied to.
    epsilon : :obj:`float`, optional
        Filter parameter at the equator, to determine cut-off where filter
        response is one-half.  That is, the scale at which the transition
        occurs.  Choose either epsilon or the filter_length_scale parameter,
        but not both.
    filter_length_scale : :obj:`float`, optional
        The filter length scale (m) at the equator.  Choose either epsilon or
        the filter_length_scale parameter, but not both.
    isotropic : :obj:`bool`, optional
        Specifically request that adjustments be made to epsilon away from the
        equator to ensure isotropy.

    Returns
    -------
        In-place operation

    Note
    ----

    References include:

    * "High-order low-pass implicit tangent filters for use in finite area
      calculations"
      WH Raymond - Monthly weather review, 1988 - journals.ametsoc.org
      (underlying formulation)
      https://doi.org/10.1175/1520-0493%281988%29116%3C2132%3AHOLPIT%3E2.0.CO%3B2
    * "A spatial filter for use in finite area calculations"
      WH Raymond, A Garder - Monthly weather review, 1988 -
      journals.ametsoc.org
      https://doi.org/10.1175/1520-0493%281988%29116%3C0209%3AASFFUI%3E2.0.CO%3B2
    * "A review of recursive and implicit filters"
      WH Raymond, A Garder - Monthly weather review, 1991 -
      journals.ametsoc.org
      (matrix inversion-cyclic case)
      https://doi.org/10.1175/1520-0493%281991%29119%3C0477%3AARORAI%3E2.0.CO%3B2
    * "High-order, high-pass implicit filters for evaluating information
      within finite areas" WH Raymond - Monthly Weather Review, 1989 -
      journals.ametsoc.org
      https://doi.org/10.1175/1520-0493%281989%29117%3C2772%3AHOHPIF%3E2.0.CO%3B2
    * "Diffusion and numerical filters" WH Raymond - Monthly weather review,
      1994 - journals.ametsoc.org
      https://doi.org/10.1175/1520-0493%281994%29122%3C0757%3ADANF%3E2.0.CO%3B2


    Warning
    -------
    - The filter is isotropic only where specified and the above criteria met.
      A warning is issued where no isotropy adjustments are made.
      Isotropic adjustment method is inherited from the CAP and will likely be
      replaced with a more generalised approach in future.

    """

    def regular_points(arr):
        uvals = np.unique(np.diff(arr))
        return ants.utils.ndarray.allclose(
            uvals[0],
            uvals,
            tolerance=ants.config.CONFIG["ants_tolerance"][
                "raymond_filter_isotropy_tolerance"
            ],
        )

    x_coord = source.coord(axis="x")
    y_coord = source.coord(axis="y")
    periodic = x_coord.circular
    src_crs = x_coord.coord_system

    do_isotropic = (
        isotropic
        and regular_points(x_coord.points)
        and regular_points(y_coord.points)
        and ants.utils.cube.is_global(source)
        and src_crs.as_ants_crs() == ants.coord_systems.UM_SPHERE.crs
    )
    if isotropic and not do_isotropic:
        msg = (
            "Currently unable to ensure isotropic filtering.  See the "
            "function docstring for conditions supported."
        )
        raise RuntimeError(msg)
    if filter_length_scale and epsilon:
        raise RuntimeError(
            'Expecting "filter_length_scale" OR "epsilon" but ' "not both."
        )

    if filter_length_scale:
        # Derive epsilon from the provided filter length scale in the NS and EW
        # direction.
        delta_lambda = x_coord.points[1] - x_coord.points[0]
        earth_radius = src_crs.semi_major_axis

        epsilon = _raymond.filter_length_scale_to_eps(
            filter_length_scale, delta_lambda, earth_radius
        )
        eps_rows = epsilon
        if do_isotropic:
            eps_rows = _raymond.filter_length_scale_to_eps_iso(
                filter_length_scale, y_coord.points, delta_lambda, earth_radius
            )
    else:
        eps_rows = epsilon
        if do_isotropic:
            eps_rows = _raymond.epsilon_iso(
                (y_coord.points.size, x_coord.points.size), epsilon
            )
    eps_columns = epsilon

    data, _ = ants.analysis._merge.horizontal_grid_reorder(source)
    _raymond.raymond_filter_ndarray(data, eps_columns, eps_rows, periodic)

    # The filter does not generate bit comparable results each time so we
    # override with 32bit values as a workaround.
    data[:] = data.astype("float32")
