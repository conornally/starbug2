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
import glob
import os

from astropy.io.fits import PrimaryHDU, HDUList
from astropy.table import Table

from starbug2.core.main_components.artificial_stars import compile_results
from starbug2.constants import ExitStates, HELP_STRINGS, Modes
from starbug2.core.star_bug_config import StarBugMainConfig
from starbug2.jwst_support.initialise_psf_data import (
    init_starbug_for_jwst, generate_psf)
from starbug2.utilities import param
from starbug2.utilities.misc import generate_runscript
from starbug2.utilities.utils import (
    printf, get_version, usage, p_error, get_data_path, split_file_name,
    export_region, combine_tables, fill_nan)


def starbug_one_time_runs(config: StarBugMainConfig) -> ExitStates:
    """
    executes any one time runs as required by the config
    :param config: the main config
    :return: the exit states of any executed
    :rtype: ExitStates
    """
    if config.show_version:
        printf(get_version())

    if config.show_help:
        usage(__doc__, verbose=config.verbose_logs)

        if config.do_star_detection:
            p_error(HELP_STRINGS[Modes.DETECTION])
        if config.do_bgd_estimate:
            p_error(HELP_STRINGS[Modes.BACKGROUND])
        if config.do_aperture_photometry:
            p_error(HELP_STRINGS[Modes.APP_HOT])
        if config.do_photometry_routine:
            p_error(HELP_STRINGS[Modes.PSFP_HOT])
        if config.do_matching:
            p_error(HELP_STRINGS[Modes.MATCH_OUTPUTS])
        return ExitStates.EXIT_EARLY

    # Load parameter files for onetime runs
    if not config.update_param:
        parameter_file: str | None
        if (parameter_file := config.param_file) is None:
            if os.path.exists("./starbug.param"):
                parameter_file = "starbug.param"
            else:
                parameter_file = None

        config.load_params(parameter_file)
    else:
        param.update_param_file(config.param_file)
        return ExitStates.EXIT_SUCCESS

    output: int | float | str
    if _output := config.output_file:
        _output: int | float | str
        output = _output
    else:
        output = '.'

    # One time run commands

    # Initialise or update starbug
    if config.execute_jwst_initialisation:
        init_starbug_for_jwst(config)

    # Generate a single PSF
    if config.generate_psf:
        if config.got_valid_psf_generation_params():
            filter_string: str | None = config.custom_filter
            assert filter_string is not None
            detector: str | None = config.detector_name
            psf_size: int = config.psf_fit_size
            if psf_size is not None:
                printf(
                    "Generating PSF: %s %s (%d)\n" %
                    (filter_string, detector, psf_size))
            else:
                printf(
                    "Generating PSF: %s %s\n" %
                    (filter_string, detector))
            psf: PrimaryHDU | None = generate_psf(
                filter_string, detector=detector, fov_pixels=psf_size)
            if psf:
                name: str = (
                    "%s%s.fits" %
                    (filter_string, "" if detector is None else detector))
                printf("--> %s\n" % name)
                d_name: str = get_data_path()
                psf.writeto(os.path.join(d_name, name), overwrite=True)
            else:
                p_error("PSF Generation failed :(\n")
        else:
            # noinspection SpellCheckingInspection
            p_error(
                "Unable to generate PSF. Set filter with '-s FILTER=FXXX and "
                "Set detector name with '-s DET_NAME=XXX and "
                "Set psf_fit_size with '-s PSF_SIZE=XXX'\n")

    # Generate a run script
    if config.generate_run:
        generate_runscript(config.fits_images, "starbug2 ")
        if not config.fits_images:
            p_error("no files included to create runscript with\n")

    # Generate a region from a table
    if config.generate_region:
        file_name: str | None = config.region_file
        if file_name and os.path.exists(file_name):
            table: Table = Table.read(file_name, format="fits")
            _, name, _ = split_file_name(file_name)
            name: str
            export_region(
                table, colour=config.region_colour,
                scale_radius=config.region_scale,
                region_radius=config.region_radius,
                x_col=config.region_x_column_name,
                y_col=config.region_y_column_name,
                wcs=config.region_uses_wcs,
                f_name="%s/%s.reg" % (output, name))
            printf("generating region --> %s/%s.reg\n" % (output, name))

    # generate local param file as requested
    if config.generate_local_param_file:
        config.do_generate_local_param_file()

    return ExitStates.EXIT_SUCCESS


def ast_one_time_runs(config: StarBugMainConfig) -> ExitStates:
    """
    executes any artificial star tess one time runs.
    :param config: the main Starbug config
    :return: exit state of any runs
    :rtype: ExitStates
    """
    if config.show_ast_help:
        usage(__doc__, verbose=config.verbose_logs)
        return ExitStates.EXIT_EARLY

    if config.ast_recover:
        f_names: list[str] | None
        if not config.fits_images:
            # noinspection SpellCheckingInspection
            f_names = glob.glob("sbast-autosave*.tmp")
        else:
            f_names = [a for a in config.fits_images if os.path.exists(a)]
        if f_names:
            printf("Recovery Mode:\n-> %s\n" % ("\n-> ".join(f_names)))
            raw: Table | None = Table()
            for f_name in f_names:
                f_name: str
                read_table: Table | None = Table.read(f_name)
                if read_table is None:
                    p_error(f"failed to read table at path {f_name}")
                    return ExitStates.EXIT_FAIL
                raw = combine_tables(raw, read_table)
            results: dict[str, HDUList]
            assert raw is not None
            if (results := compile_results(
                    fill_nan(raw), config, plot_ast="recovered.pdf")):
                printf("-> successful recovery!\n--> %s\n" % (
                    f_name := "recovered.fits"))
                for result_file_name in results.keys():
                    results[result_file_name].writeto(
                        f"{result_file_name}{f_name}", overwrite=True)
            else:
                p_error("something went wrong\n")
        else:
            p_error("No files found to recover\n")
    return ExitStates.EXIT_SUCCESS
