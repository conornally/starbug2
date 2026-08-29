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

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
import starbug2
from starbug2.constants import ExitStates, TableColumn, HeaderTags
from starbug2.core.plot import load_style, plot_test, plot_inspect_source
from starbug2.core.star_bug_config import StarBugMainConfig
from starbug2.utilities.utils import p_error, warn, parse_cmd, usage
from astropy.io.fits import PrimaryHDU, ImageHDU, BinTableHDU

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
        config.generate_plot_get_opt_definitions())

    _, argv = parse_cmd(argv)
    config.populate_params(
        argv, short_definition, long_definition, config.PLOT_FLAG_MAP)
    return config


def plot_one_time_runs(config: StarBugMainConfig) -> ExitStates:
    """
    Handles initialisation routines such as style sheet distribution or
     help menu warnings.
    """
    if config.show_plot_help:
        usage(__doc__, verbose=bool(config.verbose_logs))
        if config.inspect_parameter:
            p_error(str(fn_pinspect.__doc__))
        return ExitStates.EXIT_EARLY

    # Only throw an error if files are missing when they are explicitly
    # required
    if not config.test_mode and len(config.fits_images) == 0:
        p_error(
            "Error: Image or catalogue argument targets must be provided.\n")
        return ExitStates.EXIT_EARLY

    if config.plot_style is not None:
        load_style(str(config.plot_style))

    if config.dark_mode:
        load_style(f"{starbug2.__path__[0]}/extras/dark.style")

    return ExitStates.EXIT_SUCCESS


def fn_pinspect(
        inspect_string: str,
        images: list[fits.PrimaryHDU | fits.ImageHDU |
                     fits.BinTableHDU | None] | None = None,
        tables: list[Table] | None = None) -> plt.Figure | None:
    """
    Plot cutouts at a source position across a range of images.

    This requires a source list to be loaded, a list of image
    files, and the source catalogue number to be given. This will
    take the form::

        $~ starbug2-plot -I CN123 source_list.fits image*.fits

    :param inspect_string: the string to inspect
    :type inspect_string: str
    :param images: The list of FITS image HDUs to cut out from
    :type images: list
    :param tables: The source list containing coordinates matching a
                   'Catalogue_Number'
    :type tables: list of astropy.Table
    :return: The output figure object containing rendered cutouts
    :rtype: matplotlib.pyplot.Figure or None
    """
    fig: plt.Figure | None = None

    if inspect_string and images and tables and len(tables) > 0:
        if (TableColumn.CAT_NUM in tables[0].colnames
                and inspect_string in tables[0][TableColumn.CAT_NUM]):
            i: np.ndarray = np.where(
                tables[0][TableColumn.CAT_NUM] == inspect_string)[0]
            fig = plot_inspect_source(tables[0][i], images)
    else:
        p_error(
            f"Must include the source {TableColumn.CAT_NUM}, "
            f"a list of images and a source list \n"
        )
    return fig


def plot_main(argv: list[str]) -> ExitStates | None:
    """
    Main runtime entry path configuration structure for
    data visualisation parsing loops.
    """
    warn("Still in development\n\n")

    config: StarBugMainConfig = starbug_parse_argv(argv)

    load_style(f"{starbug2.__path__[0]}/extras/starbug.style")

    if plot_one_time_runs(config) == ExitStates.EXIT_EARLY:
        return ExitStates.EXIT_EARLY

    images: list[PrimaryHDU | ImageHDU | BinTableHDU | None] = []
    tables: list[Table] = []

    for arg in config.fits_images:
        if os.path.exists(arg):
            fp: fits.HDUList = fits.open(arg)
            _filter: str = fp[0].header.get(HeaderTags.FILTER)

            # Use type tracking alias explicitly during extraction loop blocks
            hdu: PrimaryHDU | ImageHDU | BinTableHDU | None = None
            for hdu in fp:
                if hdu is None:
                    continue

                if hdu.header.get(HeaderTags.EXT) == HeaderTags.IMAGE:
                    images.append(hdu)
                    break
                if hdu.header.get(HeaderTags.EXT) == HeaderTags.BIN_TABLE:
                    tables.append(Table(hdu.data))
                    break
            if hdu is not None:
                hdu.header[HeaderTags.FILTER] = _filter

    fig: plt.Figure | None = None

    if config.test_mode:
        ax: plt.Axes
        fig, ax = plt.subplots(1, figsize=(3, 2.5))
        plot_test(ax)

    inspect_string: str | None = config.inspect_parameter
    if inspect_string is not None:
        fig = fn_pinspect(inspect_string, images=images, tables=tables)

    if fig is not None:
        fig.tight_layout()
        if output := config.output_file:
            fig.savefig(str(output), dpi=300)
        else:
            plt.show()

    return ExitStates.EXIT_SUCCESS


def plot_main_entry() -> ExitStates | None:
    """Command Line package gateway binary endpoint entry pointer mapper."""
    return plot_main(sys.argv)
