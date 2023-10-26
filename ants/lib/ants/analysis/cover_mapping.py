# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import json
import warnings

import ants
import ants.fileformats.cover_mapping
import ants.utils
import iris
import numpy as np

try:
    import numba

    _NUMBA_IMPORT_ERROR = False
    optimise = numba.jit(nopython=True, nogil=True)
except Exception as _NUMBA_IMPORT_ERROR:
    numba = None
    msg = " {}\nProceeding without capabilities provided by numba."
    warnings.warn(msg.format(str(_NUMBA_IMPORT_ERROR)))

    def optimise(arg):
        # No-op - we don't have the numba library available to apply
        # optimisation.
        return arg


@optimise
def _apply_transform(barray, crosswalk, flag_values):
    # If numba is not available, raise error as performance is badly
    # effected. Cannot pass exception in, as unsupported in numba.
    if numba is None:
        raise ModuleNotFoundError("No module named 'numba'")

    assert flag_values.shape[0] == crosswalk.shape[0]

    rarray = np.zeros(
        (crosswalk.shape[1], barray.shape[0], barray.shape[1]), dtype=np.int8
    )
    # Loop over the data
    for j in range(barray.shape[0]):
        for i in range(barray.shape[1]):
            # Determine which src class
            for k in range(crosswalk.shape[0]):
                # src contributes?
                if barray[j, i] == flag_values[k]:
                    # Contributes to which target classes?
                    for l in range(crosswalk.shape[1]):
                        contrib = crosswalk[k, l]
                        if contrib > 0:
                            rarray[l, j, i] = rarray[l, j, i] + contrib
    return rarray


def fetch_lct_slices(source, um_tile_ids):
    """
    Fetch the slices corresponding to the specified JULE tile ID.

    Given an iris cube, derive a set of slices for this cube corresponding to
    this JULE tile ID.  That is, indexing the pseudo_level mapped dimension.

    Parameters
    ----------
    source : :class:`~iris.cube.Cube`
        Source cube which has land cover types as defined by a pseudo level
        coordinate.
    um_tile_ids : one or more int
        JULES tile IDs.  This corresponds to the pseudo_level value(s) desired.

    Returns
    -------
    : tuple
        A tuple containing slice objects.

    Example usage::

        c3_grass = 3
        c3_grass_slices = fetch_veg_slice(cube, c3_grass)
        c3_grass.data[c3_grass] = ...

    """
    um_tile_ids_iter = um_tile_ids
    if not hasattr(um_tile_ids, "__iter__"):
        um_tile_ids_iter = [um_tile_ids]
    pseudo_level = source.coord("pseudo_level")
    pseudo_level_points = pseudo_level.points.tolist()
    index = [pseudo_level_points.index(um_tile_id) for um_tile_id in um_tile_ids_iter]

    # iris requires tuples...
    index = tuple(index)
    if not hasattr(um_tile_ids, "__iter__"):
        # If we supply a single index we don't want to return an iterable slice
        # as that means that the object being sliced will have an extra
        # dimension which the user would not expect.
        index = index[0]

    slices = [slice(None)] * source.ndim
    ps_dim = source.coord_dims(pseudo_level)
    if len(ps_dim) != 1:
        msg = (
            "Expecting 1D pseudo level coordinate describing JULES "
            "classes, got {}D coord."
        ).format(len(ps_dim))
        raise RuntimeError(msg)
    slices[source.coord_dims(pseudo_level)[0]] = index
    return tuple(slices)


def normalise_fractions(source):
    """
    Normalisation of fractions, ensuring the sum of the fractions is 1.

    Normalisation effectively works by filling missing data by dividing equally
    amongst all types where data > 0 such that the ratios between fractions
    remain the same.  Where there are no class fractions at a given point, it
    will remain 0.  Fractions outside the range [0, 1] are considered anomalous
    and pulled into that range before normalisation occurs.

    Parameters
    ----------
    source : `iris.cube.Cube`
        Source land cover type fraction, with pseudo-level coordinate
        representing the classes.

    Warning
    -------
    Mask is not altered by this function.

    """
    # Ensure fractions are within a suitable range
    # Ignore numpy invalid value warnings due to presence of nans' which we
    # expect after a regrid.
    with np.errstate(invalid="ignore"):
        source.data[source.data < 0] = 0
        source.data[source.data > 1] = 1

    pseudo_level = source.coord("pseudo_level")
    if pseudo_level.ndim != 1:
        msg = "Expecting a 1D pseudo_level coordinate not {}D."
        raise RuntimeError(msg.format(pseudo_level.ndim))

    pdim = source.coord_dims(pseudo_level)[0]

    transpose_indx = [pdim] + [x for x in range(source.ndim) if x != pdim]
    data = np.asarray(ants.utils.ndarray.transposed_view(source.data, transpose_indx))

    with np.errstate(invalid="ignore"):
        non_zero = (data > 0).sum(axis=0) > 0
    if not non_zero.all():
        warnings.warn(
            "Locations present with no classification fraction, "
            "ignoring such locations."
        )
    non_zero_data = data[..., non_zero]
    source_sum = non_zero_data.sum(axis=0)
    adjustments = 1.0 - source_sum
    data[..., non_zero] = (
        adjustments / (1 - adjustments) * non_zero_data
    ) + non_zero_data


class CoverMapper(object):
    """
    Generic transformation class for applying a cover map table, transforming
    between source and target classifications.

    """

    def __init__(self, source_types, target_types, cover_map):
        """
        Returns an object describing the mapping between one set of surface
        types with another.

        Parameters
        ----------
        source_types : iterable
            An array-like object describing the source types (mapping the first
            dimension of cover_map).
        target_types : iterable
            An array-like object describing the target types (mapping the
            second dimension of cover_map).  These should be the pseudo level
            values of the surface tiles.
        cover_map : iterable
            An array-like object which describes the mapping between 'n' source
            types to 'm' target types.  Values are converted to integer
            percentages, where fractions are lost.

        """
        cover_map = np.array(cover_map)
        source_types = self._init_source_types(source_types)
        target_types = np.array(target_types, dtype="int")

        if cover_map.shape != (source_types.size, target_types.size):
            msg = (
                "The cover map has {} source mappings to {} target "
                "classifications while there are {source} expected source "
                "types and {target} expected target types."
            )
            raise RuntimeError(
                msg.format(
                    *cover_map.shape, source=source_types.size, target=target_types.size
                )
            )

        # Check cover_map totals and convert to integer percentages
        cover_map_totals = cover_map.sum(1)
        if np.allclose(cover_map_totals, 1):
            # Fractional definition
            cover_map = np.rint((cover_map * 100)).astype("int8")
        elif np.allclose(cover_map_totals, 100):
            # Percentage definition
            cover_map = np.rint(cover_map).astype("int8", copy=False)
        else:
            msg = (
                "Incorrectly defined cover map table, contributions do "
                "not add up to 100%"
            )
            raise RuntimeError(msg)

        # Adjust percentages to ensure that they add up to 100 by adjusting
        # largest fraction.
        adjustments = 100 - cover_map.sum(1)
        for cc, adjustment in zip(cover_map, adjustments):
            cc[cc.argmax()] += adjustment

        self._source_types = source_types
        self._target_types = target_types
        self._cover_map = cover_map

    def _init_source_types(self, source_types):
        return np.char.lower(np.array(source_types))

    @property
    def source_types(self):
        """Return an :class:`numpy.ndarray` describing the source types."""
        return self._source_types

    @property
    def target_types(self):
        """Return an :class:`numpy.ndarray` describing the target types."""
        return self._target_types

    @property
    def cover_map(self):
        """
        Return an :class:`numpy.ndarray` describing the cover map.

        An array-like object which describes the mapping between 'n' source
        types to 'm' target types.

        """
        return self._cover_map

    def _get_reordered_cover_map(self, source_types):
        """Return re-ordered a cover_map based on the source_types provided."""
        source_types = self._init_source_types(source_types)
        nsource_types = self._source_types
        cover_map = self._cover_map

        if source_types.size != nsource_types.size:
            raise ValueError("Underspecified source types provided.")
        missing = set(nsource_types).symmetric_difference(set(source_types))
        if missing:
            raise ValueError(
                "Source types do not match cover table "
                "description: {}".format(missing)
            )

        # Ensure that the cover map is ordered as the source
        if (source_types != nsource_types).any():
            # Source flag values do not match the order of the cover table
            # order.  Re-order the cover table.
            order = [np.argmax(nsource_types == meaning) for meaning in source_types]
            nsource_types = nsource_types[order]
            cover_map = cover_map[order]
        return nsource_types, cover_map

    def reorder_cover_map(self, source_types):
        """
        Re-order the cover_map in-place based on the source_types provided.

        Parameters
        ----------
        source_types : iterable
            An array-like object describing the source types (mapping the first
            dimension of cover_map).

        """
        nsource_types, cover_map = self._get_reordered_cover_map(source_types)
        self._source_types = nsource_types
        self._cover_map = cover_map

    def apply_transform(self, data, flag_values):
        """
        Apply transform to the provided classification array.

        Classification types are assumed to vary along axis 0
        Any masked values are set to flag 0.

        Parameters
        ----------
        data : :class:`np.ndarray`
            Data consisting of classifications.  Expecting 8-bit integers as
            its datatype.
        flag_values : :class:`np.ndarray`
            1D integer array representing the classes of 'data'. These values
            are assumed to follow the ordering of the source types.  If not,
            re-order the flag_values, or re-order the cover_map via
            :meth:`~CoverMapper.reorder_cover_map`.

        Returns
        -------
        : :class:`np.ndarray`
            Resulting data mapped onto the target classifications.

        Note
        ----
        Transform applied using 8-bit integer calculations.
        Any unmapped types are considered missing data values and so the mask
        is dropped.  Performing a normalisation of values is then suggested
        to fill these missing data values, effectively a nearest neighbour
        fill.

        """
        cover_map = self._cover_map

        if np.can_cast(flag_values, np.int8, "same_kind"):
            flag_values = flag_values.astype(np.int8, copy=False)
        if np.can_cast(cover_map, np.int8, "same_kind"):
            cover_map = cover_map.astype(np.int8, copy=False)

        data_array = data
        if np.ma.isMaskedArray(data):
            data_array = data.data
        res = _apply_transform(data_array, cover_map, flag_values)

        # Apply mask where appropriate.
        if np.ma.isMaskedArray(data) and np.ma.is_masked(data):
            res = np.ma.array(res, copy=False)
            res.mask = np.repeat(data.mask[np.newaxis], res.shape[0], axis=0)
        return res

    @classmethod
    def load(cls, filename):
        """Load the transform from a json file on disk."""
        with open(filename, "r") as fh:
            trans_dic = json.load(fh)
        try:
            transform = cls(
                trans_dic["source"], trans_dic["target"], trans_dic["cover_map"]
            )
        except KeyError as err:
            msg = (
                'Valid covermap files should include "source", "target" '
                'and "cover_map"'
            )
            err.args = ["{}. {}".format(list(err.args)[0], msg)]
            raise
        return transform

    def save(self, filename):
        """Save the transform to disk in the form of a json."""
        cover_map = self.cover_map.tolist()
        source_types = self.source_types.tolist()
        target_types = self.target_types.tolist()

        result = {
            "cover_map": cover_map,
            "source": source_types,
            "target": target_types,
        }
        with open(filename, "w") as outfile:
            json.dump(result, outfile)


class SCTTransformer(object):
    def __init__(self, transform):
        """
        Surface cover type classification transformer.

        Return an object capable of transforming between source and target
        surface classification types.

        Parameters
        ----------
        transform : :class:`~ants.analysis.cover_mapping.CoverMapper`
            Transform object which describes the mapping between the source and
            target surface classification types.

        """
        self._transform = transform

    def _group_whole_contributions(self, source, flag_values, flag_meanings):
        # Optimisation by grouping - removal of redundant types.
        transform = self._transform
        rm_list = []
        for ind in range(len(transform.target_types)):
            one_loc = np.where(transform.cover_map[:, ind] == 100)[0]
            if len(one_loc) > 1:
                ref_loc = one_loc[0]
                redundant_loc = one_loc[1:]
                rm_list.extend(redundant_loc)
                for flg_rm in redundant_loc:
                    source.data[source.data == flag_values[flg_rm]] = flag_values[
                        ref_loc
                    ]

        if rm_list:
            accept_list = set(range(len(transform.source_types))).difference(
                set(rm_list)
            )
            accept_list = list(accept_list)
            flag_values = flag_values[accept_list]
            flag_meanings = flag_meanings[accept_list]
            cover_map = transform.cover_map[accept_list, :]
            ants.fileformats.cover_mapping.set_flag_arrays(
                source, flag_values, flag_meanings
            )
            transform = CoverMapper(flag_meanings, transform.target_types, cover_map)
        return source, transform

    def _make_cube(self, dat, source):
        target_cube = iris.cube.Cube(dat)
        ants.utils.cube.inherit_metadata(target_cube, source)
        sx, sy = ants.utils.cube.horizontal_grid(source, dim_coords=True)
        target_cube.add_dim_coord(sx, source.coord_dims(sx)[0] + 1)
        target_cube.add_dim_coord(sy, source.coord_dims(sy)[0] + 1)

        points = np.array(self._transform.target_types)
        sc = iris.coords.AuxCoord(points, long_name="pseudo_level")
        target_cube.add_aux_coord(sc, 0)
        return target_cube

    def __call__(self, source, target_grid=None):
        """
        Transform the surface classification types of source to the target
        types onto the target grid.

        Parameters
        ----------
        source : :class:`~iris.cube.Cube`
            Source cube which has surface cover type data (defined by flag
            values).
        target_grid : :class:`~iris.cube.Cube`, optional
            Target grid definition.  If not provided, no regridding will take
            place.

        Returns
        -------
        : :class:`~iris.cube.Cube`
            Source data, transformed to the surface classification types
            described by the transform, on the target grid specified.

        """
        # Retrieve transform metadata from cube
        flag_values, flag_meanings = ants.fileformats.cover_mapping.get_flag_arrays(
            source
        )

        # Check and re-order transform based on cube metadata
        self._transform.reorder_cover_map(flag_meanings)

        # Optimise by taking advantage of full fraction contributions.
        source, optim_transform = self._group_whole_contributions(
            source, flag_values, flag_meanings
        )

        # Apply transform
        flag_values, flag_meanings = ants.fileformats.cover_mapping.get_flag_arrays(
            source
        )
        data = optim_transform.apply_transform(source.data, flag_values)
        source = self._make_cube(data, source)
        if target_grid:
            source = ants.analysis.mean(source, target_grid)
        source.data = source.data / 100.0

        # Normalise fractions to effectively fill missing data and account for
        # precision loss in type casting.
        normalise_fractions(source)

        # Convert to int8 where it makes sense in order to save on storage.
        if (
            target_grid is None
            and (np.unique(optim_transform.cover_map) == [0, 100]).all()
        ):
            source.data = source.data.astype("int8")
        return source
