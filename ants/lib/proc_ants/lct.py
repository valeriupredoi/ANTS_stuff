# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants
import iris
import numpy as np
from ants.analysis.cover_mapping import fetch_lct_slices, normalise_fractions


def _get_transposed_source(source):
    # Return a transposed array view so that we can utilise broadcasting rules
    # with masks.
    pseudo_level = source.coord("pseudo_level")
    pdim = source.coord_dims(pseudo_level)[0]
    transpose_indx = list(range(source.ndim))
    transpose_indx.pop(pdim)
    transpose_indx += [pdim]

    return ants.utils.ndarray.transposed_view(source.data, transpose_indx)


def set_whole_fraction_ice(source):
    """
    Set ice fractions to [0, 1].

    Parameters
    ----------
    source : `iris.cube.Cube`
        Source land cover type fraction, with pseudo-level coordinate
        representing the classes.

    See Also
    --------
    :func:`~ants.analysis.cover_mapping.normalise_fractions` : for details of
        how fractions are normalised.

    """
    # Ensure that all fractions add up to 1 before applying the threshold.
    # This function currently assumed this to be true as the crosswalk
    # calculation already normalises before returning a result.
    # normalise_fractions(source)
    threshold = 0.5

    pseudo_level = source.coord("pseudo_level")
    if pseudo_level.ndim != 1:
        msg = 'Expecting 1D "peudo_level" coordinate, got {}D'
        raise ValueError(msg.format(pseudo_level.ndim))

    ice_tile_id = 9
    slices = fetch_lct_slices(source, ice_tile_id)
    ice_data = source.data[slices]

    data = np.asarray(_get_transposed_source(source))

    # Set data above threshold.
    mask = ice_data >= threshold
    data[mask] = 0
    ice_data[mask] = 1

    # Set data below threshold - this isn't the inverse of the mask above as
    # we may have 'nan' values present which we don't want to lose.
    mask = ice_data < threshold
    ice_data[mask] = 0
    normalise_fractions(source)


def remove_ocean_level(source, min_frac=0.5):
    """
    Remove ocean level from the provided source.

    Removing the ocean level (pseudo_level == 0) occurs by redistributing cells
    with less than 50% ocean to the land fractions present to ensure they add
    to 100%

    Parameters
    ----------
    source : :class:`iris.cube.Cube`
        Source land cover type fraction, with pseudo-level coordinate
        representing the classes.
    min_frac : :obj:`float`, optional
        Land fraction below which land is removed and all vegetation fractions
        are ignored. Set this to 0.0 if supplying your own land sea mask.

    Returns
    -------
    : :class:`iris.cube.Cube`
        Vegetation fractions with the ocean level removed.

    """
    land_constraint = iris.Constraint(pseudo_level=lambda cell: cell > 0)

    # Redistribute cells which are < 0.5 area ocean coverage to the other
    # classifications.
    ocean_tile_id = 0
    slices = fetch_lct_slices(source, ocean_tile_id)
    ocean_data = source.data[slices]

    mask = ocean_data < (1.0 - min_frac)
    ocean_data[mask] = 0
    normalise_fractions(source)
    source = source.extract(land_constraint)

    # Set ocean (>= ocean) as masked.
    if not np.ma.isMaskedArray(source.data):
        source.data = np.ma.asarray(source.data)

    # This isn't the inverse of the mask above as we may have 'nan' values
    # present which we don't want to lose.
    mask = ocean_data >= (1.0 - min_frac)
    source.data[:, mask] = np.ma.masked

    # Extract lsm cube from this, inherit unknown values from the lct in the
    # form of a mask (lct mask becomes the lsm data and lct nan becomes the lcm
    # mask).
    lsm = source.slices_over("pseudo_level").next()
    lsm = lsm.copy(
        np.ma.array(
            (~np.ma.getmaskarray(lsm.data)).astype("int8"),
            mask=np.isnan(np.asarray(lsm.data)),
        )
    )
    lsm.remove_coord("pseudo_level")
    lsm.attributes["valid_min"] = 0
    lsm.attributes["valid_max"] = 1
    lsm.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi("m01s00i030")
    lsm.rename("land_binary_mask")
    return source, lsm


def remove_non_glacial_ice(lct):
    """
    Replace isolated ice and replace with soil.

    - Identify ice seed points by those locations with 5x5 surrounding ice
      cells.
    - Find contiguous ice regions using these seed points.
    - Replace ice with soil for those not matching these contiguous regions.

    """
    ice_id = 9
    ice_level = ants.analysis.cover_mapping.fetch_lct_slices(lct, ice_id)
    soil_id = 8
    soil_level = ants.analysis.cover_mapping.fetch_lct_slices(lct, soil_id)
    ice_level_cube = lct[ice_level]

    # Find ice seed points (5x5 areas of ice)
    mask = np.ones((5, 5))
    ice_seed = ants.analysis.MooreNeighbourhood(
        ice_level_cube, window_mask=mask
    ).all_equal_value(1)

    # Use seed points (above) to find contiguous ice regions and construct a
    # mask which excludes these contiguous ice regions.
    ice_data = lct.data[ice_level]
    replace_ice_mask = ice_data > 0  # all ice assigned for removal.
    while ice_seed.any():
        y, x = np.unravel_index(ice_seed.argmax(), ice_seed.shape)
        ice_contiguous = ants.analysis.find_similar_region(
            ice_data,
            (y, x),
            extended_neighbourhood=True,
            wraparound=lct.coord(axis="x").circular,
        )
        # Remove contiguous region from seed array.
        ice_seed[ice_contiguous] = False
        # Remove contiguous ice from ice removal mask.
        replace_ice_mask[ice_contiguous] = False

    # Replace ice with soil.
    soil_data = lct.data[soil_level]
    soil_data[replace_ice_mask] = 1
    ice_data[replace_ice_mask] = 0
