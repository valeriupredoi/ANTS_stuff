#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
"""
Pre-processor for the CCI (ESA) source for land cover types
***********************************************************

The following steps are taken:

1. Correct/populate CCI metadata:

  - Set the flag_values and flag_meanings.
  - Correct bounds on both 'y' and 'x' axes.

2. Merge the raw IGBP dataset over the Antarctic region (this region is
   invalid in the CCI dataset as the ice does not extend beyond the coast).
3. Fill missing data using the spiral search.
4. Set Aral sea and lake Victoria as ocean (rather than in-land water).
   The motivation to be more closely aligned to the ocean landsea mask.

"""
import logging
import warnings

import ants
import ants.io.save as save
import ants.utils.cube
import proc_ants.lct_preproc_cci as lct_preproc_cci
from proc_ants.lakes import fill_lakes

_LOGGER = logging.getLogger(__name__)


def fill_nemo_mask_lakes(cci_cube, nemo_flag):
    """
    Identify lakes as ocean as per nemo mask.

    Ocean type includes: Open ocean, Caspian Sea, Black Sea, Aral Sea,
    Lake Victoria and Great Lakes (Lake Superior, Michigan, Huron, Ontario,
    Erie).  Note that Caspian sea and Black sea is already defined as
    "sea_ocean_water", so remains unaltered by this processing.

    """
    lakes = [
        ("victoria", (-1.1, 32.9)),
        ("south aral sea", (45.0, 58.5)),
        ("superior", (47.7, -87.5)),
        ("huron", (44.78, -82.21)),
        ("michigan", (43.86, -87.09)),
        ("ontario", (43.85, -77.77)),
        ("erie", (42.25, -81.16)),
        ("north aral sea", (46.3, 61.0)),
    ]
    for lake, seed in lakes:
        try:
            cci_cube = fill_lakes(
                cci_cube, lake, seed, "water_bodies", nemo_flag, constrain=False
            )
        except ants.exceptions.FloodfillError:
            warnings.warn(f"Lake {lake} already has type {nemo_flag}")
    return cci_cube


def fill_ostia_lakes(cci_cube, ostia_flag):
    """
    Identify water bodies in the CCI that are to be set as "resolved lake" type

    We identify the largest lakes in OSTIA as per ARC-Lake database, filling
    them as a new classification type ("resolved lake" or -47).  All other
    water bodies remain as they are.  Top 70 lakes (order by area minimum) from
    the ARC-Lake database are included.  Also included are those lakes with no
    minimum area size defined.  Those chosen have instead a polygon area greater
    than the area min of the smallest lake of above, yathkyed (1449km2).

    """
    # Lakes/reservoirs ordered by decreasing area minimum km2 (high->low) or
    # decreasing area polygon where there is no area minimum available.
    lakes = [
        ("tanganyika", (-6.07, 29.46)),  # area_min:32000.00km2; area_poly:32820.50km2
        ("baikal", (53.63, 108.14)),  # area_min:31500.00km2; area_poly:31924.60km2
        ("bear", (65.91, -121.30)),  # area_min:31326.00km2; area_poly:30530.10km2
        ("slave", (62.09, -114.37)),  # area_min:28568.00km2; area_poly:27816.30km2
        ("winnipeg", (52.12, -97.25)),  # area_min:24387.00km2; area_poly:23809.30km2
        ("malawi", (-11.96, 34.59)),  # area_min:22490.00km2; area_poly:29251.50km2
        ("ladoga", (60.84, 31.39)),  # area_min:18130.00km2; area_poly:17539.10km2
        ("balkhash", (45.91, 73.95)),  # area_min:17000.00km2; area_poly:17458.80km2
        ("onega", (61.90, 35.35)),  # area_min:9700.00km2; area_poly:9608.10km2
        ("nicaragua", (11.57, -85.36)),  # area_min:8150.00km2; area_poly:7851.50km2
        ("titicaca", (-15.92, -69.30)),  # area_min:8030.00km2; area_poly:8270.60km2
        ("athabasca", (59.10, -109.96)),  # area_min:7935.00km2; area_poly:7781.60km2
        ("reindeer", (57.19, -102.27)),  # area_min:6650.00km2; area_poly:5596.60km2
        ("smallwood", (54.19, -64.31)),  # area_min:6475.00km2; area_poly:5610.40km2
        ("turkana", (3.53, 36.08)),  # area_min:6400.00km2; area_poly:7785.40km2
        ("issyk-kul", (42.46, 77.25)),  # area_min:6240.00km2; area_poly:6258.90km2
        ("albert", (1.67, 30.91)),  # area_min:5590.00km2; area_poly:5401.90km2
        ("nern", (58.88, 13.22)),  # area_min:5580.00km2; area_poly:5550.50km2; vanern
        ("nettilling", (66.42, -70.28)),  # area_min:5542.00km2; area_poly:5064.70km2
        ("winnipegosis", (52.37, -100.05)),  # area_min:5375.00km2; area_poly:5167.10km2
        ("nipigon", (49.80, -88.55)),  # area_min:4848.00km2; area_poly:4476.90km2
        # ? ("manitoba", (50.99, -98.80)),  # area_min:4625.00km2; area_poly:4790.90km2
        ("qinghai", (36.89, 100.18)),  # area_min:4460.00km2; area_poly:4449.70km2
        # ? ("kyoga", (1.50, 33.01)),  # area_min:4430.00km2; area_poly:1727.70km2
        # ? ("kwania", (1.72, 32.65)),  # area_min:4430.00km2; area_poly:561.90km2
        ("mweru", (-9.01, 28.74)),  # area_min:4350.00km2; area_poly:4944.80km2
        # ? ("woods", (49.38, -94.91)),  # area_min:4350.00km2; area_poly:4167.70km2
        # ? ("peipus", (58.41, 27.59)),  # area_min:4300.00km2; area_poly:3552.80km2
        # ? ("taymyr", (74.42, 101.76)),  # area_min:4000.00km2; area_poly:4200.10km2
        # ? ("khanka", (44.94, 132.42)),  # area_min:4000.00km2; area_poly:4088.10km2
        # ? ("dubawnt", (63.13, -101.44)),  # area_min:3833.00km2; area_poly:3628.50km2
        # ? (
        # ?     "savern",
        # ?     (38.66, 42.98),
        # ? ),  # area_min:3740.00km2; area_poly:3537.90km2;
        # # We are actually asking for lake "van" but it's geometry is associated
        # # with "savern" under Natural Earth...
        # ("tana", (11.95, 37.31)),  # area_min:3600.00km2; area_poly:3034.70km2
        # ("uvs", (50.33, 92.81)),  # area_min:3350.00km2; area_poly:3421.50km2
        ("amadjuak", (64.99, -71.13)),  # area_min:3115.00km2; area_poly:3033.70km2
        # ("mirim", (-32.89, -53.25)),  # area_min:2970.00km2; area_poly:3958.90km2
        # ("hungtze", (33.34, 118.53)),  # area_min:2700.00km2; area_poly:1524.90km2
        # ("wollaston", (58.30, -103.33)),  # area_min:2681.00km2; area_poly:2272.00km2
        # ("alakol", (46.11, 81.75)),  # area_min:2650.00km2; area_poly:2802.10km2
        # ("hovsgol", (51.02, 100.48)),  # area_min:2620.00km2; area_poly:2741.40km2
        # ("iliamna", (59.56, -154.90)),  # area_min:2590.00km2; area_poly:2623.20km2
        # ("nam", (30.71, 90.66)),  # area_min:2500.00km2; area_poly:1933.60km2
        # ("mistassini", (50.82, -73.81)), # area_min:2335.00km2; area_poly:2162.70km2
        # ("nueltin", (60.25, -99.40)),  # area_min:2279.00km2; area_poly:2011.90km2
        # ("buenos aires", (-46.66, -72.50)),
        # # area_min:2240.00km2; area_poly:1848.20km2
        # ("kivu", (-2.04, 29.23)),  # area_min:2220.00km2; area_poly:2370.70km2
        # ("tai", (31.21, 120.24)),  # area_min:2210.00km2; area_poly:2398.50km2
        # ("edward", (-0.39, 29.61)),  # area_min:2150.00km2; area_poly:2252.50km2
        # ("vattern", (58.33, 14.57)),  # area_min:1910.00km2; area_poly:1847.20km2
        # ("baker", (64.13, -95.28)),  # area_min:1887.00km2; area_poly:1780.30km2
        # ("ziling", (31.77, 88.95)),  # area_min:1860.00km2; area_poly:1640.90km2
        # ("chiquita", (-30.74, -62.61)),  # area_min:1850.00km2; area_poly:2011.90km2
        # ("okeechobee", (26.95, -80.86)),  # area_min:1810.00km2; area_poly:1436.80km2
        # ("martre", (63.33, -117.91)),  # area_min:1776.00km2; area_poly:1684.00km2
        # ("har us", (48.06, 92.30)),  # area_min:1760.00km2; area_poly:1193.50km2
        # ("puruvesi", (61.77, 29.02)),  # area_min:1760.00km2; area_poly:1086.60km2
        # ("orivesi", (62.35, 29.59)),  # area_min:1760.00km2; area_poly:958.70km2
        # ("haukivesi", (62.10, 28.52)),  # area_min:1760.00km2; area_poly:830.60km2
        # ("eskimo", (69.10, -132.76)),  # area_min:1676.00km2; area_poly:1146.10km2
        # ("hulun", (48.97, 117.38)),  # area_min:1590.00km2; area_poly:2166.50km2
        # (
        #     "tengiz",
        #     (50.44, 68.90),
        #     False,
        # ),  # area_min:1590.00km2; area_poly:1382.60km2;
        # # tengiz fetches the SELETYTENIZ lake geometry so we force a non match.
        # ("yathkyed", (62.69, -98.07)),  # area_min:1449.00km2; area_poly:1320.40km2
        # ("itaparica", (-10.18, -42.01)),  # area_min:--km2; area_poly:8779.10km2
        # ("volta", (7.63, 0.11)),  # area_min:--km2; area_poly:7490.70km2
        # ("great salt lake", (41.20, -112.50)),  # area_min:--km2; area_poly:5965.80km2
        # ("kuybyshevskoye", (54.54, 48.65)),  # area_min:--km2; area_poly:5003.30km2
        # ("urmia", (37.64, 45.49)),  # area_min:--km2; area_poly:4963.40km2
        # ("kariba", (-17.23, 27.60)),  # area_min:--km2; area_poly:4958.50km2
        # ("bratsk", (55.93, 101.86)),  # area_min:--km2; area_poly:4476.50km2
        # ("zaysan", (48.70, 83.44)),  # area_min:--km2; area_poly:4470.50km2
        # ("cabora  bassa", (-15.73, 31.63)),  # area_min:--km2; area_poly:4362.60km2
        # ("rybinkskoye", (58.49, 38.13)),  # area_min:--km2; area_poly:3926.60km2
        # ("nasser", (22.86, 32.58)),  # area_min:--km2; area_poly:3817.90km2
        # ("tucurui", (-4.57, -49.49)),  # area_min:--km2; area_poly:3485.20km2
        # ("grande", (53.70, -74.87), False),  # area_min:--km2; area_poly:3107.60km2
        # ("volgogradskoye", (50.35, 45.85)),  # area_min:--km2; area_poly:2697.90km2
        # ("grande", (53.86, -76.73), False),  # area_min:--km2; area_poly:2572.60km2
        # ("boeng tonle chhma", (12.81, 104.15)),
        # # area_min:--km2; area_poly:2569.90km2
        # ("tsimlyanskoye", (48.05, 42.98)),  # area_min:--km2; area_poly:2242.40km2
        # ("southern indian", (57.14, -98.61)),  # area_min:--km2; area_poly:2227.40km2
        # ("tapajos", (-2.88, -55.14)),  # area_min:--km2; area_poly:2210.50km2
        # (
        #     "buhayrat ath tharthar",
        #     (34.11, 43.17),
        # ),  # area_min:--km2; area_poly:2119.40km2
        # ("poyang", (29.25, 116.06)),  # area_min:--km2; area_poly:2108.10km2
        # ("kakhovskoye", (47.27, 33.95)),  # area_min:--km2; area_poly:2083.90km2
        # ("vilyuyskoye", (62.73, 111.16)),  # area_min:--km2; area_poly:2076.70km2
        # ("bangweulu", (-11.19, 29.76)),  # area_min:--km2; area_poly:2049.10km2
        # ("manicouagan", (51.35, -69.13)),  # area_min:--km2; area_poly:2034.70km2
        # ("syvash", (45.96, 34.74)),  # area_min:--km2; area_poly:1967.60km2
        # ("rukwa", (-7.84, 32.16)),  # area_min:--km2; area_poly:1965.80km2
        # ("kremenshugskoye", (49.28, 32.62)),  # area_min:--km2; area_poly:1964.90km2
        # ("zeyskoye", (54.26, 127.80)),  # area_min:--km2; area_poly:1964.00km2
        # ("ijsselmeer", (52.66, 5.42)),  # area_min:--km2; area_poly:1962.70km2
        # ("mai-ndombe", (-2.14, 18.32)),  # area_min:--km2; area_poly:1954.90km2
        # ("chad", (13.08, 14.52)),  # area_min:--km2; area_poly:19346.60km2
        # ("chany", (54.83, 77.39)),  # area_min:--km2; area_poly:1918.00km2
        # ("krasnoyarskoye", (54.84, 90.94)),  # area_min:--km2; area_poly:1893.20km2
        # ("argyle", (-16.43, 128.77)),  # area_min:--km2; area_poly:1749.00km2
        # ("ust-ilimskoye", (57.15, 102.32)),  # area_min:--km2; area_poly:1657.70km2
        # ("williston", (55.95, -123.91)),  # area_min:--km2; area_poly:1654.00km2
        # ("kamskoye", (58.80, 56.26)),  # area_min:--km2; area_poly:1648.20km2
        # ("blommesteinmeer", (4.72, -55.07)),  # area_min:--km2; area_poly:1631.40km2
        # ("seul", (50.40, -92.05)),  # area_min:--km2; area_poly:1610.80km2
        # ("kossou", (7.55, -5.68)),  # area_min:--km2; area_poly:1566.00km2
        # ("lakes sakakawea", (47.81, -102.32)),  # area_min:--km2; area_poly:1492.00km2
        # ("saimmaa", (61.39, 28.20)),  # area_min:--km2; area_poly:1466.10km2
    ]

    for lake in lakes:
        try:
            lake_name, seed, use_geom = lake
        except ValueError:
            lake_name, seed = lake
            use_geom = True
        try:
            cci_cube = fill_lakes(
                cci_cube,
                lake_name,
                seed,
                ["water_bodies", "sea_ocean_water"],
                ostia_flag,
                constrain=True,
                use_lake_geom=use_geom,
            )
        except ants.exceptions.FloodfillError:
            warnings.warn(f"Lake {lake_name} already has type {ostia_flag}")
    return cci_cube


def load_data(cci_filepath, igbp_filepath):
    """
    Load data relevant to pre-processing the ESA CCI data source.

    Parameters
    ----------
    cci_filepath : str
        CCI source data filepath.
    igbp_filepath : str
        Filepath for the IGBP source (for filling missing data in the CCI).

    Returns
    -------
    : tuple(:class:`~iris.cube.Cube`, :class:`~iris.cube.Cube`)
        Representing CCI source and IGBP source respectively.

    """
    cci_cube = ants.load_cube(cci_filepath, "land_cover_lccs")
    ants.utils.cube.set_crs(cci_cube)

    igbp_cube = ants.load_cube(igbp_filepath)

    return cci_cube, igbp_cube


def main(
    source_filepath,
    igbp_filepath,
    output_filepath,
    nemo_flag,
    ostia_flag,
    use_new_saver,
):
    """
    Pre-processing ESA CCI top level call function.

    Perform the neccesary pre-processing to produce a CCI source suitable for
    deriving JULES land cover types.

    Parameters
    ----------
    source_filepath : str
        CCI source data filepath.
    glacial_geom_filepath : str
        Non-polar glacial geometries filepath.
    igbp_filepath : str
        Filepath for the IGBP source (for filling missing data in the CCI,
        specifically Antarctica).
    output_filepath : str
        Output file path for the resulting ancillary.

    Returns
    -------
    : :class:`~iris.cube.Cube`
        Pre-processed ESA CCI source.

    """
    cci_cube, igbp_cube = load_data(source_filepath, igbp_filepath)
    _LOGGER.info("load completed")
    lct_preproc_cci.update_cci_metadata(cci_cube)
    _LOGGER.info("metadata update completed")
    cci_cube = lct_preproc_cci.merge_igbp(cci_cube, igbp_cube)
    _LOGGER.info("merge with igbp completed")

    cci_cube.data
    _LOGGER.info("data touch completed")
    filler = ants.analysis.FillMissingPoints(cci_cube)
    filler(cci_cube)
    _LOGGER.info("spiral completed")

    cci_cube = fill_nemo_mask_lakes(cci_cube, nemo_flag)
    _LOGGER.info("Marking lakes as ocean completed")

    cci_cube = fill_ostia_lakes(cci_cube, ostia_flag)
    _LOGGER.info("Re-classification of lakes as 'resolved lakes' completed")

    if use_new_saver:
        save.netcdf(cci_cube, output_filepath)
    else:
        ants.save(cci_cube, output_filepath)

    return cci_cube


def _get_parser():
    parser = ants.AntsArgParser()
    parser.add_argument(
        "--igbp-source",
        type=str,
        help="IGBP source filepath which is used to fill " "Antarctica in the CCI.",
        required=True,
    )
    parser.add_argument(
        "--nemo-flag",
        type=str,
        help="Flag to assign for NEMO lakes.",
        default="sea_ocean_water",
    )
    parser.add_argument(
        "--ostia-flag",
        type=str,
        help="Flag to assign for OSTIA lakes.",
        default="resolved_lake",
    )
    return parser


if __name__ == "__main__":
    args = _get_parser().parse_args()
    main(
        args.sources,
        args.igbp_source,
        args.output,
        args.nemo_flag,
        args.ostia_flag,
        args.use_new_saver,
    )
