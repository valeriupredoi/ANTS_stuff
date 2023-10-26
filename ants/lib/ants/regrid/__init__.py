# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
"""
ANTS regridding provides capability which extend beyond what is currently
provided by iris for the convenience of ancillary generation.
Therefore, the user is referred to :mod:`iris.analysis` for regridding
capability provided by iris.  ANTS provides:

* :mod:`ants.regrid.rectilinear` a set of rectilinear horizontal
  regridding/interpolation approaches.
* :mod:`ants.regrid.interpolation` a set of vertical points-based approaches.
* :mod:`ants.regrid.esmf` regridding schemes for ESMF.

The reader is referred to the module documentation for further details.
For further details see the user guide.

"""
import logging
import sys

from ants.config import CONFIG

from . import _ugrid, esmf, interpolation, rectilinear

_LOGGER = logging.getLogger(__name__)


class GeneralRegridder(object):
    def __init__(
        self, src_grid, target_grid, horizontal_regridder=None, vertical_regridder=None
    ):
        """
        General regridder abstracting away horizontal and vertical regridding.

        Parameters
        ----------
        src_grid : :class:`~iris.cube.Cube`
           Defining the source grid.
        target_grid : :class:`~iris.cube.Cube`
            Defining the target grid.
        horizontal_regridder : :obj:`str`, optional
            Horizontal regridder callable.
        vertical_regridder : :obj:`str`, optional
            Vertical regridder callable.

        """
        self._horizontal_regridder = horizontal_regridder
        self._vertical_regridder = vertical_regridder
        if horizontal_regridder is None and vertical_regridder is None:
            raise AttributeError(
                "At least one of horizontal or vertical re-grid schemes must be "
                "provided."
            )

    def __call__(self, cube):
        """
        Regrid both vertical and horizontally where appropriate.

        Parameters
        ----------
        cube : :class:`~iris.cube.Cube`
            Source to be regridded.

        Returns
        -------
        : :class:`~iris.cube.Cube`
            Redridded result.

        """
        res = cube
        if self._vertical_regridder:
            res = self._vertical_regridder(res)
        if self._horizontal_regridder:
            res = self._horizontal_regridder(res)

        return res


class GeneralRegridScheme(object):
    """
    Abstract away the concept of horizontal and vertical regridding by
    providing a general scheme that handles both under the hood.

    """

    def __init__(self, horizontal_scheme=None, vertical_scheme=None):
        """
        General scheme which handles both vertical and horizontal regrid.

        The GeneralRegridScheme is useful to define a regrid method(s) and
        allow this regridding to be overridden after the fact via a
        configuration file where necessary.  In the case where a fixed
        regridding scheme is wanted and no override is to be allowed, please
        use the regridding scheme directly.

        Parameters
        ----------
        horizontal_scheme : :obj:`str`, optional
            Name of horizontal regridding scheme to use.
            Default regridding scheme is None.
        vertical_scheme : :obj:`str`, optional
            Name of vertical regridding scheme to use.
            Default regridding scheme is None.

        """
        regridder_scheme = (
            CONFIG["ants_regridding_horizontal"]["scheme"] or horizontal_scheme
        )
        regridder = regridder_scheme
        if isinstance(regridder_scheme, str):
            regridder_scheme = (
                getattr(rectilinear, regridder, None)
                or getattr(sys.modules["iris.analysis"], regridder, None)
                or getattr(esmf, regridder, None)
                or getattr(_ugrid, regridder)
            )
            regridder = regridder_scheme()
        self._horizontal_scheme = regridder

        regridder_scheme = (
            CONFIG["ants_regridding_vertical"]["scheme"] or vertical_scheme
        )
        extrapolation_mode = (
            CONFIG["ants_regridding_vertical"]["extrapolation_mode"] or None
        )
        regridder = regridder_scheme
        if isinstance(regridder_scheme, str):
            regridder_scheme = getattr(interpolation, regridder, None)
            if extrapolation_mode is not None:
                regridder = regridder_scheme(extrapolation=extrapolation_mode.lower())
            else:
                regridder = regridder_scheme()
        self._vertical_scheme = regridder
        _LOGGER.info(repr(self))

    def regridder(self, src_grid, target_grid):
        """
        Creates a GeneralRegridder to regrid from the source to target grid.

        Parameters
        ----------
        src_grid : :class:`~iris.cube.Cube`
           Defining the source grid.
        target_grid : :class:`~iris.cube.Cube`
            Defining the target grid.

        Returns
        -------
        : callable
           Callable with the interface `callable(cube)`

           where `cube` is a cube with the same grid as `src_grid`
           that is to be regridded to the `target_grid`.

        """
        hregridder = None
        if self._horizontal_scheme is not None:
            hregridder = self._horizontal_scheme.regridder(src_grid, target_grid)
        vregridder = None
        if self._vertical_scheme is not None:
            vregridder = self._vertical_scheme.regridder(src_grid, target_grid)
        return GeneralRegridder(src_grid, target_grid, hregridder, vregridder)

    def __repr__(self):
        return "{}(horizontal_scheme={!r}, vertical_scheme={!r})".format(
            self.__class__.__name__, self._horizontal_scheme, self._vertical_scheme
        )
