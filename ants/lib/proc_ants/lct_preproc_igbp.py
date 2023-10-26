#!/usr/bin/env python2.7
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants
from ants.fileformats.cover_mapping import set_flag_arrays


def _missing_metadata(cube):
    """
    Fill missing metadata in the IGBP cube.

    The following corrections are made:
    - flag_values.
    - flag_meanings.
    - x and y bounds.
    - Populate missing coordinate system (WGS84).
    - Provide a suitable long_name (IGBP land classification).

    Parameters
    ----------
    cube : :class:`~iris.cube.Cube`
        IGBP source cube.

    Returns
    -------
    : :class:`~iris.cube.Cube`
        IGBP source with corrections and with populated metadata.

    """
    cube.long_name = "IGBP land classification"
    cube.var_name = None
    cube.attributes["source"] = (
        "Data derived from USGS IGBP land " "classification data"
    )
    # see https://lta.cr.usgs.gov/glcc/globdoc2_0#app2
    flag_meanings = (
        "ocean evergreen_needleleaf_forest evergreen_broadleaf_forest "
        "deciduous_needleleaf_forest deciduous_broadleaf_forest "
        "mixed_forest closed_shrublands open_shrublands woody_savannas "
        "savannas grasslands permanent_wetlands croplands urban_and_built-up "
        "cropland_or_natural_vegetation_mosaic snow_and_ice "
        "barren_or_sparsely_vegetated water_bodies"
    )
    set_flag_arrays(cube, range(18), flag_meanings)
    cube.units = 1
    ants.utils.cube.set_crs(cube, ants.coord_systems.WGS84_GEODETIC)
    # Correct grid definition problems.
    y_coord = cube.coord(axis="y")
    bounds = y_coord.bounds.copy()
    bounds[0, 0], bounds[-1, -1] = -90, 90
    y_coord.bounds = bounds

    x_coord = cube.coord(axis="x")
    bounds = x_coord.bounds.copy()
    bounds[0, 0], bounds[-1, -1] = -180, 180
    x_coord.bounds = bounds


def pre_process(igbp, bats=None):
    """
    Pre-process the IGBP.

    Pre-processing the igbp includes differentiating between ocean and in-land
    water by using the 'bats' dataset and secondly, metadata for the IGBP is
    corrected.

    Parameters
    ----------
    igbp : :class:`~iris.cube.Cube`
        IGBP source cube.
    bats : :class:`~iris.cube.Cube`, optional
        BATS source cube.  If not provided, only metadata correction of the
        IGBP will occur.

    """

    def operation(igbp, bats):
        ocean = 15
        igbp.data[bats.data == ocean] = 0
        return igbp

    if bats is not None:
        operation(igbp, bats)
    _missing_metadata(igbp)
    return igbp
