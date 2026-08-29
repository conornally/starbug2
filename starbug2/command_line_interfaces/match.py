"""Copyright (C) 2026 UKATC

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>."""
import os
import sys
from typing import Any

import numpy as np
from astropy.table import Table, vstack
from astropy.units import Quantity
from starbug2.utilities import utils
from starbug2.constants import (
     STAR_BUG_MIRI, NIRCAM, MATCH_COLS, TableColumn, ExitStates, HeaderTags)
from starbug2.utilities.filters import STAR_BUG_FILTERS
from starbug2.matching.band_match import BandMatch
from starbug2.matching.cascade_match import CascadeMatch
from starbug2.matching.exact_value_match import ExactValueMatch
from starbug2.matching.generic_match import GenericMatch
from starbug2.utilities.misc import parse_mask
from starbug2.core.star_bug_config import StarBugMainConfig
from starbug2.utilities.utils import parse_cmd, usage
import photutils

# Force photutils to strictly return standard QTables globally
photutils.future_column_names = True


def starbug_parse_argv(argv: list[str]) -> StarBugMainConfig:
    """
    Organise the sys argv line into options, values and arguments

    :param argv: the arguments
    :return: the config class
    :rtype: StarBugMainConfig
    """
    config: StarBugMainConfig = StarBugMainConfig()
    short_definition: str
    long_definition: list[str]
    short_definition, long_definition = (
        config.generate_match_get_opt_definitions())

    _, argv = parse_cmd(argv)
    config.populate_params(
        argv, short_definition, long_definition, config.MATCH_FLAG_MAP)
    return config


def match_full_band_match(
        tables: list[Table],
        d_threshold: np.ndarray, bridge_col: str) -> Table:
    """
    Handles fallback deprecated band-matching configurations across diverse
     detectors.
    """
    utils.p_error("THIS NEEDS A TEST\n")
    to_match: dict[int, list[Table]] = {
        NIRCAM: [],
        STAR_BUG_MIRI: []
    }
    _col_names: list[str] = [TableColumn.RA, TableColumn.DEC, TableColumn.FLAG]
    band_matcher: BandMatch = BandMatch(threshold=d_threshold)

    for tab in tables:
        tab: Table
        filter_string: str = str(tab.meta.get(HeaderTags.FILTER))
        to_match[STAR_BUG_FILTERS[filter_string].instr].append(tab)
        _col_names += [filter_string, f"e{filter_string}"]

    matched: Table
    if to_match[NIRCAM] and to_match[STAR_BUG_MIRI]:
        utils.printf("Detected NIRCam to MIRI matching\n")
        nir_cam_matched: Table = band_matcher.band_match(
            to_match[NIRCAM], col_names=_col_names)
        miri_matched: Table = band_matcher.band_match(
            to_match[STAR_BUG_MIRI], col_names=_col_names)

        # noinspection SpellCheckingInspection
        load: utils.Loading = utils.Loading(
            len(miri_matched),
            msg=f"Combining NIRCAM-MIRI({np.array2string(d_threshold)}g\")"
        )

        mask: np.ndarray
        if bridge_col:
            mask = np.isnan(nir_cam_matched[bridge_col])
            utils.printf("-> bridging catalogues with %s\n" % bridge_col)
        else:
            mask = np.full(len(nir_cam_matched), False)

        m: GenericMatch = GenericMatch(threshold=d_threshold, load=load)
        full: Any = m((nir_cam_matched[~mask], miri_matched))
        matched = m.finish_matching(full)
        matched.remove_column(TableColumn.NUM)
        matched = vstack((matched, nir_cam_matched[mask]))
    else:
        matched = band_matcher.band_match(tables, col_names=_col_names)

    return matched


def match_main(argv: list[str]) -> ExitStates:
    """
    Main runtime processing loop for executing cross-catalogue astronomical
    source coordinate matching.
    """
    config: StarBugMainConfig = starbug_parse_argv(argv)

    if config.show_match_help:
        usage(__doc__, verbose=config.verbose_logs) # noqa
        return ExitStates.EXIT_SUCCESS

    p_file: str | None = config.param_file
    if not p_file:
        config.param_file = (
            "./starbug.param" if os.path.exists("./starbug.param") else None)

    tables: list[Table] = []
    for f_name in config.fits_table:
        t: Table | None = utils.import_table(f_name, verbose=True)
        if t is not None:
            tables.append(t)

    masks: list[np.ndarray | None] = []
    if raw := config.mask_eval:
        masks = [parse_mask(raw, t) for t in tables]
        for m in masks:
            try:
                print(str(m), sum(m), len(m)) # noqa
            except (TypeError, NameError, ImportError):
                print(m) # noqa

    if len(tables) > 1:
        col_names: list[str] = list(MATCH_COLS)

        if config.extra_match_columns is not None:
            col_names += [
                name for name in config.extra_match_columns.split()
                if name not in col_names
            ]
        d_threshold: Quantity = config.match_threshold_arc_sec_as_an_arc_sec
        error_column: str = config.error_col

        average_table: Table
        output_table: Table | None = None
        matcher: GenericMatch

        if config.band_deprecated:
            average_table = match_full_band_match(
                tables, config.match_threshold_arc_sec_as_an_array,
                config.bridge_band_column)
            config.full_run = True
        else:
            if config.do_band_processing:
                band_threshold: np.ndarray = (
                    config.match_threshold_arc_sec_as_an_array)

                filter_string: list[str]
                if config.custom_filter != "":
                    filter_string = str(config.custom_filter).split(',')
                else:
                    filter_string = utils.remove_duplicates(
                        [utils.find_filter(t) for t in tables])

                matcher = BandMatch(
                    threshold=band_threshold, fltr=filter_string,
                    verbose=config.verbose_logs
                )
            elif config.do_cascade:
                matcher = CascadeMatch(
                    threshold=d_threshold, colnames=col_names,
                    verbose=config.verbose_logs
                )
            elif config.generic_mode:
                matcher = GenericMatch(
                    threshold=d_threshold, col_names=col_names,
                    verbose=config.verbose_logs
                )
            elif config.exact_match:
                matcher = ExactValueMatch(
                    value=TableColumn.CAT_NUM, colnames=None,
                    verbose=config.verbose_logs
                )
            else:
                matcher = GenericMatch(
                    threshold=d_threshold, verbose=config.verbose_logs
                )
                config.full_run = True

            if config.verbose_logs:
                print("\n%s" % matcher)

            output_table = matcher.match(tables, join_type="or", mask=masks)
            average_table = matcher.finish_matching(
                output_table,
                num_thresh=config.exposure_count_threshold,
                zp_mag=config.zero_point_magnitude,
                error_column=error_column
            )

        output: str | None = config.output_file
        if output is None or output == '.':
            output = utils.combine_file_names(
                [name for name in config.fits_images], n_mismatch=100)
            if output is None:
                return ExitStates.EXIT_FAIL

        d_name: str
        f_name: str
        ext: str
        d_name, f_name, ext = utils.split_file_name(output)

        suffix: str = ""
        if config.full_run and output_table is not None:
            utils.export_table(
                output_table, f_name="%s/%sfull.fits" % (d_name, f_name))
            utils.printf("-> %s/%sfull.fits\n" % (d_name, f_name))
            suffix = "match"

        if average_table:
            utils.export_table(
                average_table, "%s/%s%s.fits" % (d_name, f_name, suffix))
            utils.printf("-> %s/%s%s.fits\n" % (d_name, f_name, suffix))

        return ExitStates.EXIT_SUCCESS

    elif len(tables) == 1:
        return ExitStates.EXIT_EARLY
    else:
        utils.p_error("No tables loaded for matching.\n")
        return ExitStates.EXIT_EARLY


def match_main_entry() -> ExitStates:
    """StarbugII-match entry path map setup routing wrapper."""
    return match_main(sys.argv)
