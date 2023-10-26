# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
from ants.fileformats.cover_mapping import set_flag_arrays


def missing_metadata(cube):
    cube.long_name = "ITE land classification"
    cube.attributes["source"] = (
        "Data derived from CEH/ITE land " "classification of Great Britain 1990."
    )
    # see https://www.ceh.ac.uk/services/land-cover-map-1990
    flag_meanings = (
        "sea_or_estuary inland_water beach_and_coastal_bare "
        "saltmarsh grass_heath mown_or_grazed_turf "
        "meadow_or_verge_meadow_or_verge_or_semi-natural "
        "rough_or_marsh_grass moorland_grass open_shrub_moor "
        "dense_shrub_moor bracken dense_shrub_heath scrub_or_orchard "
        "deciduous_woodland coniferous_woodland upland_bog tilled_land "
        "ruderal_weed suburban_or_rural_development "
        "continuous_urban inland_bare_ground felled_forest lowland_bog "
        "open_shrub_heath"
    )
    set_flag_arrays(cube, range(1, 26), flag_meanings)
    cube.units = 1


def pre_process_source(cube):
    missing_metadata(cube)
    return cube
