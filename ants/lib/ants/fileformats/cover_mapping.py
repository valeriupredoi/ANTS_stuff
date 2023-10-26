# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants
import ants.analysis.cover_mapping
import numpy as np


def set_flag_arrays(cube, flag_values, flag_meanings):
    """
    Set CF compliant flag_values and flag_meanings.

    Parameters
    ----------
    cube : :class:`iris.cube.Cube`
        Source cube to specify flag values and meanings.
    flag_meanings : string or iterable
        Flag meanings where when specified by a single string, a comma, colon,
        space and semicolon are acceptable delimiters.
    flag_meanings : list or np.ndarray
        Flag values.

    """
    if ~np.equal(np.trunc(flag_values), flag_values).any():
        raise ValueError("Flag values cannot be represented as integers")
    flag_values = np.asarray(flag_values, int).tolist()

    if isinstance(flag_meanings, str):
        delimiter = [",", ";", ":"]
        dfound = None
        for delim in delimiter:
            if delim in flag_meanings:
                if not dfound:
                    dfound = delim
                else:
                    msg = (
                        "More than one potential delimiter found in " '"flag_meanings"'
                    )
                    raise ValueError(msg)
        flag_meanings = [cont for cont in flag_meanings.split(dfound or " ") if cont]
    flag_meanings = [meaning.strip() for meaning in flag_meanings]

    if len(flag_meanings) != len(flag_values):
        msg = "Flag values and flag meanings are not the same length"
        raise ValueError(msg)

    flag_meanings = " ".join(flag_meanings)
    cube.attributes["flag_values"] = flag_values
    cube.attributes["flag_meanings"] = flag_meanings


def get_flag_arrays(cube):
    """
    Get flag_values and flag_meanings from the provided cube.

    Fetch flag_values and flag_meanings from a cube in a standardised manner
    including a consistent datatype and formatting.

    Parameters
    ----------
    cube : :class:`iris.cube.Cube`
        Source cube to fetch flag values and meanings.

    Returns
    -------
    : tuple(numpy.ndarray, numpy.ndarray)
        The tuple represented (flag_values, flag_meanings).

    Raises
    ------
    RuntimeError
        If the number of flag_values and flag_meanings in the cube do not
        match.

    """
    flag_values = np.array(cube.attributes["flag_values"], dtype="int")
    flag_meanings = np.array(cube.attributes["flag_meanings"].lower().split())
    if flag_values.size != flag_meanings.size:
        raise RuntimeError("Missing flag value/meaning pair")
    return flag_values, flag_meanings


def load_cover_mapper(transform_path):
    """
    Load a transform object for mapping source and target land cover types.

    A transform object is returned which uses the mapping information defined
    from the specified json file.  This object then allows the mapping of
    source land cover types to the target land cover types.

    Parameters
    ----------
    transform_path : str
        Path to the json transform file.

    Returns
    -------
    : :class:`~ants.analysis.cover_mapping.CoverMapper`

    """
    return ants.analysis.cover_mapping.CoverMapper.load(transform_path)
