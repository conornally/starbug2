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
import warnings

from astropy.io.fits.verify import VerifyWarning
from astropy.io.fits.header import Header
from astropy.table import Table
from astropy.utils.exceptions import AstropyDeprecationWarning, AstropyWarning
import photutils

from starbug2.core.main_components.multi_treading_execution import (
    execute_multi_core_main, execute_one_core_run_main)
from starbug2.core.main_components.one_time_runs import starbug_one_time_runs
from starbug2.constants import (
    ExitStates, LOGO, READ_THE_DOCS_URL)
from starbug2.matching.generic_match import GenericMatch
from starbug2.core.star_bug_config import StarBugMainConfig
from starbug2.core.starbug_main import StarbugBase
from starbug2.utilities.utils import (
    combine_file_names, export_table, parse_cmd, p_error, printf, puts,
    split_file_name)

# Target-silence only the specific Photutils/Astropy deprecation noise
# without masking generic Runtime math errors globally.
warnings.filterwarnings(
    "ignore", category=AstropyDeprecationWarning)
warnings.filterwarnings(
    "ignore", message=".*contains deprecated section.*",
    category=AstropyWarning)

# Handle RuntimeWarnings elegantly: Ignore expected ones (like NaN comparisons
# during clipping), but let actual mathematical issues surface.
warnings.filterwarnings(
    "ignore", message=".*invalid value encountered.*", category=RuntimeWarning)
warnings.filterwarnings(
    "ignore", message=".*divide by zero.*", category=RuntimeWarning)

# --- FITS IO FORMATTING NOISE ---
# These suppress warnings about FITS header compliance
# (e.g., truncated comments) that do not affect scientific output.
warnings.filterwarnings(
    "ignore",
    category=VerifyWarning,
    message=".*Card is too long.*"
)

# Force photutils to strictly return standard QTables globally
photutils.future_column_names = True

# noinspection SpellCheckingInspection
"""StarbugII - JWST PSF photometry
usage: starbug2 [-ABDfGhMPSv] [-b bgdfile] [-d apfile] [-n ncores] [-o ouput]
                [-p file.param] [-s opt=val] image.fits ...
   -A  --apphot          : run aperture photometry on a source list
   -B  --background      : run background estimation
   -D  --detect          : run source detection
   -G  --geom            : calculate geometric stats on source list
   -M  --match           : match outputs from all input image files
   -P  --psf             : run psf photometry
   -S  --subbgd          : subtract background from image

   -b  --bgdfile         : load background (-bgd.fits) file
   -d  --apfile  ap.fits : load a source detection (-ap.fits) file to skip the
                           source detection step
   -f  --find            : attempt to find associated -ap -bgd files
   -h  --help            : display uasage information
   -n  --ncores      num : number of CPU cores to split process between
   -o  --output      dir : output directory
   -p  --param   a.param : load parameter file
   -s  --set      option : set value in parameter file at runtime (-s SIGSKY=3)
   -v  --verbose         : display verbose outputs

   --> Single run commands
       --init                     : Initialise Starbug (post install)
       --local-param              : Make a local copy of the default
                                    parameter file
       --update-param             : Update an out-of-date local parameter file
       --generate-psf             : Generate a single PSF. Set FILTER,
                                    DET_NAME, PSF_SIZE with -s
       --generate-region   a.fits : Make a ds9 region file with a detection
                                    file
       --generate-run      *.fits : Generate a simple run script
       --version                  : Print starbug2 version

   --> typical runs
      //Source detect on image with a parameter file
      $~ starbug2 -vD -p file.param image.fits

      //Source detect and match outputs of a list of images
      $~ starbug2 -vDM -n4 images*.fits

      //PSF photometry on an image with a source file (image-ap.fits)
      $~ starbug2 -vd image-ap.fits -BP image.fits

To see more detailed information on an option, run [OPTION] --help:
    $~ starbug2 -D --help

See https://starbug2.readthedocs.io for full documentation.

"""

# noinspection SpellCheckingInspection
sys.stdout.write("\x1b[1mlaunching \x1b[36mstarbug\x1b[0m\n")


# noinspection SpellCheckingInspection
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
        config.generate_main_get_opt_definitions())

    _, argv = parse_cmd(argv)
    config.populate_params(
        argv, short_definition, long_definition, config.MAIN_FLAG_MAP)

    # shut down loading psf if we're in one time runs
    if config.use_main_one_time_runs():
        config.unfreeze()
        config.ast_load_psf = False
        config.freeze()
    return config


def starbug_match_outputs(
        starbugs: list[StarbugBase | None], config: StarBugMainConfig) -> None:
    """
    Matching output catalogues

    :param starbugs: star bug instances
    :type starbugs: list of starbugBase or None
    :param config: the config object
    :type config: StarBugMainConfig
    :return: None
    """
    if config.verbose_logs:
        printf("Matching outputs\n")

    f_name: str | None

    # filter out any Nones.
    valid_bugs: list[StarbugBase] = [sb for sb in starbugs if sb is not None]

    # get file name
    if f_name := combine_file_names([sb.f_name for sb in valid_bugs]):
        _, name, _ = split_file_name(os.path.basename(f_name))
        name: str
        f_name = "%s/%s" % (valid_bugs[0].out_dir, name)
    else:
        f_name = "out"

    header: Header = valid_bugs[0].header()

    match: GenericMatch = GenericMatch(
        threshold=config.match_threshold_arc_sec_as_an_arc_sec,
        col_names=None,
        p_file=config.param_file)

    if config.do_star_detection or config.do_aperture_photometry:
        full: Table = match(
            [sb.detections for sb in valid_bugs], join_type="or")
        av: Table = match.finish_matching(
            full, num_thresh=config.exposure_count_threshold,
            zp_mag=config.zero_point_magnitude)

        printf("-> %s-ap*...\n" % f_name)

        # noinspection SpellCheckingInspection
        export_table(full, f_name="%s-apfull.fits" % f_name, header=header)

        # noinspection SpellCheckingInspection
        export_table(av, f_name="%s-apmatch.fits" % f_name, header=header)

    if config.do_photometry_routine:
        full: Table = match(
            [sb.psf_catalogue for sb in valid_bugs], join_type="or")
        av: Table = match.finish_matching(
            full, num_thresh=config.exposure_count_threshold,
            zp_mag=config.zero_point_magnitude)

        printf("-> %s-psf*...\n" % f_name)

        # noinspection SpellCheckingInspection
        export_table(full, f_name="%s-psffull.fits" % f_name, header=header)

        # noinspection SpellCheckingInspection
        export_table(av, f_name="%s-psfmatch.fits" % f_name, header=header)


def starbug_main(argv: list[str]) -> ExitStates:
    """
    Command-line execution orchestrator for processing astronomical image
    datasets.

    :param argv: System arguments mapping configurations and input filenames
    :type argv: list of str
    :return: System operational termination exit code status matrix
    :rtype: ExitStates
    """
    config: StarBugMainConfig = starbug_main_entry_parse(argv)
    return starbug_internal_main(config)


def starbug_internal_main(config: StarBugMainConfig) -> ExitStates:
    """
   Main control for processing astronomical image datasets.

   :param config: the starbug config object
   :type config: StarBugMainConfig.
   :return: System operational termination exit code status matrix
   :rtype: ExitStates
   """
    if config.use_main_one_time_runs():
        return starbug_one_time_runs(config)

    if config.fits_images:
        # freeze the config now to avoid writers
        config.freeze()

        puts(LOGO % READ_THE_DOCS_URL)
        exit_code: ExitStates = ExitStates.EXIT_SUCCESS
        starbugs: list[StarbugBase | None]

        if ((n_cores := config.n_cores) is None
                or n_cores == 1 or len(config.fits_images) == 1):
            starbugs = execute_one_core_run_main(config)
        else:
            starbugs = execute_multi_core_main(config)

        to_remove: list[StarbugBase | None] = []
        sb: StarbugBase | None
        for n, sb in enumerate(starbugs):
            if not sb:
                p_error("FAILED: %s\n" % config.fits_images[n])
                to_remove.append(sb)
                exit_code = ExitStates.EXIT_MIXED
        for sb in to_remove:
            starbugs.remove(sb)

        if config.do_matching and len(starbugs) > 1:
            starbug_match_outputs(starbugs, config)

    else:
        p_error("fits image file must be included\n")
        exit_code = ExitStates.EXIT_FAIL

    return exit_code


def starbug_main_entry_parse(argv: list[str]) -> StarBugMainConfig:
    """Auxiliary entry parser execution wrapper."""
    return starbug_parse_argv(argv)


def starbug_main_entry() -> int:
    """
    System binary path gateway routing console script entries.
    """
    return starbug_main(sys.argv)
