# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
"""
The `ants` package provides access to the fileformats commonly used in
ancillary generation. These include those supported by iris and additional
formats such as grid namelists and raster files.

####################
Loading Target grids
####################

The ANTS library provides a grid loading capability via
:func:`ants.load_grid`.  This function will return an iris cube that can be
used to define the target grid for processing.  The cube can have a constant
data field or a land sea mask depending on the data source.

Grid loading supports the following formats:

1. Ancillary (CAP) compliant Fortran Namelist grid definition files (further
   information can be found at :func:`ants.fileformats.namelist.load_grid`).
   Use this format in applications that do not need to use the land sea mask
   of the target field as part of their processing.
2. Any fileformat that can be interpreted by iris.
   Use this format to read a grid with a land sea mask.

############
Loading data
############

ANTS uses iris to load common fileformats and therefore has access to all
formats that are supported by iris.

In most cases you should, however, use :func:`ants.load`,
:func:`ants.load_cube` etc. to load data as cubes from files.
These ANTS functions perform additional processing on the input
data.  They will attempt to:

1. Derive the global/regional status of the input data and
   set the metadata on the cube accordingly.
2. Guess bounds of latitude and longitude for any sources without
   horizontal grid bounds.
3. Remove any forecast_reference_time and forecast_period coordinates.

The ANTS functions add extra capability, depending on the fileformat, over the
base iris functions.  You should read the documentation for the ANTS module
for the specific fileformat being used, for example
:mod:`ants.fileformats.ancil` for ancillary fields file loading.

###########
Saving data
###########

Saving via ``ants.fileformats`` is deprecated as of ANTS 1.1.0 and will be
removed in a future release. The ``--use-new-saver`` argument can be used when
running applications to test the new :mod:`ants.io.save` interface.

Like loading, ANTS save adds extra file format dependent capability over
the base iris save.  Saving the cube to NetCDF uses
:func:`ants.fileformats.netcdf.cf.save`.  Saving to an F03 UM ancillary file
(using ``saver='ancil'``) uses :func:`ants.fileformats.ancil.save` and
:func:`ants.fileformats.netcdf.cf.save` to save the cube as *both* a NetCDF
and an F03 ancillary file.  NetCDF is easier to use for verification and
diagnostic purposes, while the fields file is currently needed for the UM.

Additionally, ANTS has support for 'ukca' flavoured NetCDF, chosen by
specifying ``saver='ukca'`` (see :mod:`ants.fileformats.netcdf.ukca`).

"""
import os
import sys
import types
from datetime import datetime
from functools import wraps

import ants
import iris
import iris.fileformats
import iris.io
import numpy as np
from ants.deprecations import (
    issue_save_deprecation,
    save_deprecation_message_for_docstring,
)

from . import _grid_extract
from . import _ugrid as ugrid
from . import ancil, namelist, netcdf, pp, raster
from .raster import gdal


def load_landsea_mask(filename, land_threshold=None):
    """
    Load a landsea mask from either a landsea mask file or a landfraction file.

    Parameters
    ----------
    filename : str
        Landsea mask or landsea fraction filepath.
    land_threshold : :obj:`float`, optional
        Threshold for converting the land fraction field into a landsea mask
        field.  0.5 would mean that any fraction greater than this will be
        masked.  This argument is used when loading a land fraction field.

    Returns
    -------
    : :class:`iris.cube.Cube`
        Landsea mask cube.

    """
    try:
        # Is it a landsea mask field?
        lbm = ants.load_cube(filename, "land_binary_mask")
        lbm = lbm.copy(lbm.data.astype("bool", copy=False))
    except iris.exceptions.ConstraintMismatchError:
        try:
            # Is it a land fraction field?
            land_fraction = ants.load_cube(filename, "vegetation_area_fraction")
            lbm = land_fraction.copy(land_fraction.data > land_threshold)
            lbm.rename("land_binary_mask")
        except iris.exceptions.ConstraintMismatchError:
            # It looks like we are wanting to extract a landsea mask from some
            # other field.
            cube = ants.load(filename)[0]
            y = cube.coord(axis="y")
            x = cube.coord(axis="x")
            cube = cube.slices((y, x)).next()
            lbm = cube.copy(~np.ma.getmaskarray(cube.data))
    return lbm


def _update_history_cmd(cube):
    """
    Update the cube history attribute with timestamped commandline arguments.

    Parameters
    ----------
    cube : :class:`~iris.cube.Cube`
        Cube to modify its history attribute.

    See Also
    --------
    :func:`update_history`

    """
    metadata = ants.config.CONFIG["ants_metadata"]["history"]
    date = datetime.today().replace(microsecond=0)
    cubes = ants.utils.cube.as_cubelist(cube)
    for cc in cubes:
        items = sys.argv[:]
        items[0] = os.path.basename(items[0])
        items.append("({})".format(metadata)) if metadata else None
        ants.utils.cube.update_history(cc, " ".join(items), date)
        if cc.attributes["history"] != cubes[0].attributes["history"]:
            raise RuntimeError(
                "History attributes on cubes being saved do not match:"
                f"{cc.name()}, {cc.attributes['history']}"
                f"{cubes[0].name()}, {cubes[0].attributes['history']}"
            )


# Replace the ancillary format agent with a customised (wrapped one).
iris.fileformats.FORMAT_AGENT.add_spec(
    iris.fileformats.FormatSpecification(
        "xUM Fieldsfile (FF) pre v3.1",
        iris.fileformats.MagicNumber(8),
        0x000000000000000F,
        ancil.load_cubes,
        priority=6,
        constraint_aware_handler=True,
    )
)


iris.fileformats.FORMAT_AGENT.add_spec(
    iris.fileformats.FormatSpecification(
        "xUM Fieldsfile (FF) post v5.2",
        iris.fileformats.MagicNumber(8),
        0x0000000000000014,
        ancil.load_cubes,
        priority=5,
        constraint_aware_handler=True,
    )
)


iris.fileformats.FORMAT_AGENT.add_spec(
    iris.fileformats.FormatSpecification(
        "xUM Fieldsfile (FF) ancillary",
        iris.fileformats.MagicNumber(8),
        0xFFFFFFFFFFFF8000,
        ancil.load_cubes,
        priority=4,
        constraint_aware_handler=True,
    )
)


iris.fileformats.FORMAT_AGENT.add_spec(
    iris.fileformats.FormatSpecification(
        "xUM Fieldsfile (FF) converted " "with ieee to 32 bit",
        iris.fileformats.MagicNumber(4),
        0x00000014,
        ancil.load_cubes_32bit_ieee,
        priority=4,
        constraint_aware_handler=True,
    )
)


iris.fileformats.FORMAT_AGENT.add_spec(
    iris.fileformats.FormatSpecification(
        "xUM Fieldsfile (FF) ancillary " "converted with ieee to 32 bit",
        iris.fileformats.MagicNumber(4),
        0xFFFF8000,
        ancil.load_cubes_32bit_ieee,
        priority=4,
        constraint_aware_handler=True,
    )
)


class _GdalIdentify(iris.io.format_picker.FileElement):
    """A :class:`FileElement` that queries 'gdalinfo' for the file."""

    def get_element(self, basename, file_handle):
        result = False
        if gdal is not None:
            result = True
            gdal.UseExceptions()
            try:
                gdal.Open(file_handle.name)
            except RuntimeError as err:
                if "not recognized as a supported file format" in err.args[0]:
                    result = False
        return result


iris.fileformats.FORMAT_AGENT.add_spec(
    iris.fileformats.FormatSpecification(
        "gdal",
        _GdalIdentify(),
        True,
        raster.load_cubes,
        priority=0,
        constraint_aware_handler=False,
    )
)


iris.fileformats.FORMAT_AGENT.add_spec(
    iris.fileformats.FormatSpecification(
        "UM Post Processing file (PP)",
        iris.fileformats.MagicNumber(4),
        0x00000100,
        pp.load_cubes,
        priority=6,
        constraint_aware_handler=True,
    )
)


iris.fileformats.FORMAT_AGENT.add_spec(
    iris.fileformats.FormatSpecification(
        "UM Post Processing file (PP) little-endian",
        iris.fileformats.MagicNumber(4),
        0x00010000,
        pp.load_cubes_little_endian,
        priority=4,
        constraint_aware_handler=True,
    )
)

HORIZONTAL_NAMELIST_FORMAT = iris.fileformats.FormatSpecification(
    "Namelist horizontal definition",
    iris.fileformats.LeadingLine(),
    lambda line: any([group in str(line).lower() for group in ["&horizgrid", "&grid"]]),
    namelist.load_cap_horizontal,
    priority=4,
    constraint_aware_handler=False,
)
iris.fileformats.FORMAT_AGENT.add_spec(HORIZONTAL_NAMELIST_FORMAT)

iris.fileformats.FORMAT_AGENT.add_spec(
    iris.fileformats.FormatSpecification(
        "Namelist vertical definition",
        iris.fileformats.LeadingLine(),
        lambda line: "&vertlevs" in str(line).lower(),
        namelist.load_um_vertical,
        priority=4,
        constraint_aware_handler=False,
    )
)


def load_grid(filenames, *args, **kwargs):
    """
    Load a grid definition and return an iris cube.

    The purpose of this function is to load a single grid definition, without
    the data payload, from components derived from one or more file.  That is,
    to merge these components where possible to return a single grid
    definition.  Examples include the merging of a vertical definition with
    one defining the horizontal domain.  Surface altitude (orography) is also
    supported as a merge component.  A hybrid height coordinate will ensue.
    Where source datasets have conflicting components, an exception will be
    raised.  Where more than one source dataset has identical components, a
    single grid will still be returned.

    Parameters
    ----------
    filenames : str
        Pathname of a file, or a list of such paths that contain the grid
        definition.
    *args :
        See :func:`ants.load`
    **kwargs :
        See :func:`ants.load`

    Returns
    -------
    :class:`iris.cube.Cube`
        Cube representing the grid defined with no data disk dependence and
        negligible data memory usage (using dummy data).

    Raises
    ------
    ValueError
        Where we are unable to resolve a single grid due to a conflicting
        coordinate.
    iris.exceptions.ConstraintMismatchError
        Where no cubes have been found.

    """
    if isinstance(filenames, str):
        filenames = (filenames,)

    results = iris.cube.CubeList()

    # TODO: Once UGrid loading uses ants.load interface, this block can be
    # removed.
    #
    # Load any UGrid mesh files, and collect other filenames in regular_filenames:
    regular_filenames = list()
    for filename in filenames:
        cubes = False
        try:
            cube = ugrid.load_mesh(filename)
            cubes = iris.cube.CubeList((cube,))
        except (RuntimeError, IOError):
            # IOError generated when loading variable resolution grids,
            # RuntimeError for loading other UM-space data.
            regular_filenames.append(filename)
        if cubes:
            results.extend(cubes)

    # Use normal ants load for all regular files in one go.  This is needed
    # for variable resolution namelists - need both parts (regular and
    # variable namelists) in the namelist loader at the same time.
    if regular_filenames:
        cubes = load(regular_filenames, *args, **kwargs)
        if cubes:
            results.extend(cubes)

    # Remove dependence on disc for both UGrid and regular files.
    grid = None
    if results:
        grid = _grid_extract.extract_grid(results)
    if not grid:
        raise iris.exceptions.ConstraintMismatchError("no cubes found")
    return grid


def _customised_load(func):
    def _pseudo_level_order_removal(cube):
        # Discard '_pseudo_level_order' coord on load (an implementation detail
        # we want to hide from the user).
        # Ideally this code can be removed when we have official support from
        # iris.
        cubes = ants.utils.cube.as_cubelist(cube)
        for cube in cubes:
            try:
                cube.remove_coord("_pseudo_level_order")
            except iris.exceptions.CoordinateNotFoundError:
                pass

    def _forecast_coordinates_removal(cube):
        # Discard 'forecast_period' and 'forecast_reference_time' coords on
        # load.
        cubes = ants.utils.cube.as_cubelist(cube)
        coords = ("forecast_period", "forecast_reference_time")
        for cube in cubes:
            for coord in coords:
                try:
                    cube.remove_coord(coord)
                except iris.exceptions.CoordinateNotFoundError:
                    pass

    # Wrap the iris load function with our own set of behaviours.
    @wraps(func)
    def load_function(*args, **kwargs):
        # Ensure that we leave appropriate calling to the underlying iris load
        # function.
        cubes = func(*args, **kwargs)
        if cubes is not None:
            try:
                ants.utils.cube.derive_circular_status(cubes)
            except iris.exceptions.CoordinateNotFoundError:
                pass
            _pseudo_level_order_removal(cubes)
            _forecast_coordinates_removal(cubes)
        return cubes

    return load_function


def _customised_save(func):
    # Wrap the iris save function with our own set of format specific
    # defaults as iris enforces that every fileformat must be append-able
    # except NetCDF which is hard-coded within iris as an exception to
    # this rule - https://github.com/SciTools/iris/issues/1728
    def save_function(*args, **kwargs):
        """
        Save one or more Cubes to file (or other writeable).

        This function is a wrapper to the iris save function and provides
        identical functionality other than the following:

        1. Support for saving ancillary files (saver = 'ancil').
        2. When saving as an ancillary, NetCDF output is also produced.
        3. Any arguments passed (for example, --fill_value) will only be used
           in the generation of the netcdf file and not the ancillary file.
        4. Any fileformat living in :mod:`ants.fileformats` is prioritised
           over those available in iris.
        5. Where the saver is unspecified and cannot be deduced from the
           filename, 'ancil' is the default saver.

        Parameters
        ----------

        source : :class:`iris.cube.Cube`, :class:`iris.cube.CubeList` or \
sequence of cubes.
            One or more cubes to be saved.  Note that ANTS uses v2.3 of iris which
            uses ``extract`` with the ``strict`` argument rather than
            ``extract_cube``.
        target : str
            When given a filename, the suitable fileformat is determined where
            possible.
        saver : :obj:`string`, optional
            Specifies the save function to use. If omitted, an attempt will be
            made to determine the format.
        kwargs :
            are passed through to the save function.

        See also
        --------

        :func:`iris.save`, :mod:`ants.fileformats`

        """
        # Issue a deprecation warning for the 'ants.save' interface.
        issue_save_deprecation("ants.save")

        # Ensure that we leave appropriate calling to the underlying iris save
        # function.
        if len(args) >= 0:
            source = args[0]
            target = args[1]
            saver_spec = kwargs.get("saver") or ants.config.CONFIG["saver"]

            # Update history attribute where applicable
            _update_history_cmd(source)

            # Set saver, determine which one from filename extension
            saver = set_saver(saver_spec, target)

        # Workaround for not having to support ancillary appending (iris
        # expects append keyword).  Changes would be required in iris to allow
        # this an exception (just like the netcdf saver), being capable of
        # saving cubelists.
        # TODO: Support the append keyword or modify iris to allow more
        # generic support for savers.
        if (
            "ants.fileformats.ancil" in saver.__module__
            or "_enforced_netcdf_on_ancil_save" in saver.__name__
            or "ants.fileformats.netcdf" in saver.__module__
        ):
            kwargs.pop("saver", None)
            return saver(*args, **kwargs)
        elif "iris.fileformats.netcdf" in saver.__module__:
            kwargs.pop("saver", None)
            return netcdf.cf.save(*args, **kwargs)
        else:
            return func(*args, **kwargs)

    save_function.__doc__ = (
        save_deprecation_message_for_docstring("ants.save") + save_function.__doc__
    )
    return save_function


def _enforced_netcdf_on_ancil_save(cubes, filename, *args, **kwargs):
    netcdf.cf.save(cubes, filename + ".nc", *args, **kwargs)
    ancil.save(cubes, filename)


def set_saver(saver_spec, target):
    saver = None
    _default_saver = "ancil"
    if saver_spec is None:
        if isinstance(target, str):
            saver = iris.io.find_saver(target)
        elif isinstance(target, types.FileType):
            saver = iris.io.find_saver(target.name)
        if saver is None:
            saver = iris.io.find_saver(_default_saver)
    else:
        saver = iris.io.find_saver(saver_spec)
        if saver is None:
            msg = "Cannot save; no saver can be found associated " 'with "{}"'
            raise ValueError(msg.format(saver_spec))
    return saver


iris.io.add_saver("ancil", _enforced_netcdf_on_ancil_save)


iris.io.add_saver("ukca", netcdf.ukca.save)


save = _customised_save(iris.io.save)


load_cube = _customised_load(iris.load_cube)


load = _customised_load(iris.load)


load_cubes = _customised_load(iris.load_cubes)


load_raw = _customised_load(iris.load_raw)
