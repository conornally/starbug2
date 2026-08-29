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
import getopt
import os
import numpy as np
from astropy import units
from astropy.table import Table
from astropy.units import Quantity
from typing import Dict, Tuple, Final, Any
from parse import parse
from pathlib import Path

from starbug2.constants import (
    SCI, DEFAULT_COLOUR, HeaderTags, AP_FILE, BGD_FILE, PSF_FILE, TableColumn,
    STAR_BUG_PARAMS, DEFAULT_PSF_FILE_NAME, PROBLEMATIC_FILTER_ID,
    PROBLEMATIC_FILTER_WARNING, STARBUG_DATA_DIR,
    DEFAULT_FULL_WIDTH_HALF_MAX, DEFAULT_MIN_MAG, DEFAULT_MAX_MAG)
from starbug2.utilities.filters import FilterStruct
from starbug2.utilities.utils import p_error, get_version, warn


class StarBugMainConfig:
    # the first characters to exclude
    EXCLUDES: Final[str] = "# \t\n //"

    # A single master map linking (Short Flag, Long Flag) to the internal
    # property for MAIN
    # Format: (short_flag, long_flag, type) -> property_name
    # None for short_flag means it only has a long version
    # noinspection SpellCheckingInspection
    MAIN_FLAG_MAP: Dict[Tuple[str | None, str, Any], str] = {
        ('A', 'apphot', bool): 'do_aperture_photometry',
        ('B', 'background', bool): 'do_bgd_estimate',
        ('D', 'detect', bool): 'do_star_detection',
        ('f', 'find', bool): 'find_file',
        ('G', 'geom', bool): 'do_source_geometry',
        ('h', 'help', bool): 'show_help',
        ('M', 'match', bool): 'do_matching',
        ('P', 'psf', bool): 'do_photometry_routine',
        ('S', 'subbgd', bool): 'do_bgd_subtraction',
        ('v', 'verbose', bool): 'verbose_logs',
        ('b', 'bgdfile', str): 'background_file',
        ('d', 'apfile', str): 'ap_file',
        ('n', 'ncores', int): 'n_cores',
        ('o', 'output', str): 'output_file',
        ('p', 'param', str): 'param_file',
        ('s', 'set', str): 'set_parameter',
        (None, 'init', bool): 'execute_jwst_initialisation',
        (None, 'generate-psf', bool): 'generate_psf',
        (None, 'local-param', bool): 'generate_local_param_file',
        (None, 'generate-region', bool): 'generate_region',
        (None, 'version', bool): 'show_version',
        (None, 'generate-run', bool): 'generate_run',
        (None, 'update-param', bool): 'update_param',
        (None, 'debug', bool): 'debug_mode',
        (None, 'dev', bool): 'dev_mode',
    }

    # A single master map linking (Short Flag, Long Flag) to the internal
    # property for AST
    # Format: (short_flag, long_flag, type) -> property_name
    # None for short_flag means it only has a long version
    # noinspection SpellCheckingInspection
    AST_FLAG_MAP: Dict[Tuple[str | None, str, Any], str] = {
        ("h", "help", bool): 'show_ast_help',
        ('v', 'verbose', bool): 'verbose_logs',
        ('n', 'ncores', int): 'n_cores',
        ('p', 'param', str): 'param_file',
        ('s', 'set', str): 'set_parameter',
        ('o', 'output', str): 'output_file',
        ('N', 'ntests', int): "artificial_star_tests_count",
        ('S', 'nstars', int): "stars_per_artificial_test",
        ('R', 'recover', bool): "ast_recover",
        (None, 'autosave', int): "ast_auto_save",
        (None, 'no-background', bool): "ast_no_background",
        (None, 'no-psfphot', bool): "ast_no_psf_phot",
        (None, 'save_added_image', bool): "save_added_image",
        (None, 'save_added_image_path', str): "save_added_image_path",
        (None, 'seed', int): "ast_seed"
    }

    # noinspection SpellCheckingInspection
    MATCH_FLAG_MAP: Dict[Tuple[str | None, str, Any], str] = {
        # Boolean Switches (No arguments)
        ('B', 'band', bool): 'do_band_processing',
        ('C', 'cascade', bool): 'do_cascade',
        ('G', 'generic', bool): 'generic_mode',
        ('X', 'exact', bool): 'exact_match',
        ('e', 'error', str): 'error_col',
        ('f', 'full', bool): 'full_run',
        ('h', 'help', bool): 'show_match_help',
        ('m', 'mask', str): 'mask_eval',
        ('o', 'output', str): 'output_file',
        ('p', 'param', str): 'param_file',
        ('s', 'set', str): 'set_parameter',
        ('v', 'verbose', bool): 'verbose_logs',
        (None, 'band-depr', bool): 'band_deprecated',
        (None, 'dither', bool): 'use_dither',
    }

    # noinspection SpellCheckingInspection
    PLOT_FLAG_MAP: Dict[Tuple[str | None, str, Any], str] = {
        # Boolean Switches (No arguments)
        ('h', 'help', bool): 'show_plot_help',
        ('v', 'verbose', bool): 'verbose_logs',
        ('X', 'test', bool): 'test_mode',
        (None, 'apfile', bool): 'ap_file',
        (None, 'dark', bool): 'dark_mode',

        # Options with Arguments (Strings)
        ('I', 'inspect', str): 'inspect_parameter',
        ('o', 'output', str): 'output_file',
        ('d', 'style', str): 'plot_style',
    }

    # Comprehensive mapping configuration linking keys to internal
    # properties and types
    # noinspection SpellCheckingInspection
    MAIN_PARAM_FILE_MAP: Dict[str, Tuple[str, type]] = {
        # GENERIC
        "VERBOSE": ("verbose_logs", bool),
        "OUTPUT": ("output_file", str),
        "HDUNAME": ("hdu_name", str),
        "FILTER": ("custom_filter", str),
        # DETECTION
        "FWHM": ("full_width_half_max", float),
        "SIGSKY": ("sigma_sky", float),
        "SIGSRC": ("sigma_source", float),
        "DOBGD2D": ("do_bgd_2d", bool),
        "DOCONVL": ("do_convolution", bool),
        "CLEANSRC": ("clean_sources", bool),
        "SHARP_LO": ("sharp_cutoff_low", float),
        "SHARP_HI": ("sharp_cutoff_high", float),
        "ROUND1_HI": ("round1_cutoff_high", float),
        "ROUND2_HI": ("round2_cutoff_high", float),
        "SMOOTH_LO": ("smooth_low", float),
        "SMOOTH_HI": ("smooth_high", float),
        "RICKER_R": ("ricker_wavelet_radius", float),
        # APERTURE PHOTOMETRY
        "APPHOT_R": ("aperture_phot_radius", float),
        "ENCENERGY": ("encircled_energy_fraction", float),
        "SKY_RIN": ("sky_annulus_inner_radius", float),
        "SKY_ROUT": ("sky_annulus_outer_radius", float),
        "APCORR_FILE": ("ap_corr_file_override", str),
        # BACKGROUND ESTIMATION
        "BGD_R": ("bgd_radius", float),
        "PROF_SCALE": ("profile_scaling_factor", float),
        "PROF_SLOPE": ("profile_slope", float),
        "BOX_SIZE": ("background_box_size", int),
        "BGD_CHECKFILE": ("bgd_check_file", str),
        # PHOTOMETRY
        "AP_FILE": ("ap_file", str),
        "BGD_FILE": ("background_file", str),
        "PSF_FILE": ("psf_file_override", str),
        "USE_WCS": ("use_wcs_values", bool),
        "ZP_MAG": ("zero_point_magnitude", float),
        "CRIT_SEP": ("critical_separation", float),
        "FORCE_POS": ("force_centroid_position", bool),
        "MAX_XYDEV": ("max_xy_deviation", str),
        "PSF_SIZE": ("psf_fit_size", int),
        "GEN_RESIDUAL": ("generate_residual_image", bool),
        # SOURCE STATS
        "CALC_CROWD": ("calculate_crowding_metric", bool),
        # CATALOGUE MATCHING
        "MATCH_THRESH": ("match_threshold_arc_sec", str),
        "MATCH_COLS": ("extra_match_columns", str),
        "NEXP_THRESH": ("exposure_count_threshold", int),
        "BRIDGE_COL": ("bridge_band_column", str),
        # ARTIFICIAL STAR TESTS
        "DO_AST": ("do_artificial_star_test", bool),
        "NTESTS": ("artificial_star_tests_count", int),
        "NSTARS": ("stars_per_artificial_test", int),
        "SUBIMAGE": ("sub_image_crop_size", int),
        "MAX_MAG": ("test_magnitude_bright_limit", int),
        "MIN_MAG": ("test_magnitude_faint_limit", int),
        "PLOTAST": ("ast_plot_filename", str),
        # MISC EXTRAS
        "REGION_COL": ("region_colour", str),
        "REGION_SCAL": ("region_scale", bool),
        "REGION_RAD": ("region_radius", int),
        "REGION_XCOL": ("region_x_column_name", str),
        "REGION_YCOL": ("region_y_column_name", str),
        "REGION_WCS": ("region_uses_wcs", bool),
        # generate psf extras
        "DET_NAME": ("detector_name", str),
        # generation of region tables
        "REGION_TAB": ("region_file", str),
        "PARAM_TAG": ("param_tag", str),

        # --- NEW PARAM FILE SHORT-CIRCUITS FOR STEPS & FLOW CONTROLS ---
        "RUN_APPHOT": ("do_aperture_photometry", bool),
        "RUN_BGD_EST": ("do_bgd_estimate", bool),
        "RUN_DETECT": ("do_star_detection", bool),
        "RUN_GEOM": ("do_source_geometry", bool),
        "RUN_MATCH": ("do_matching", bool),
        "RUN_PSFPHOT": ("do_photometry_routine", bool),
        "RUN_BGDSUB": ("do_bgd_subtraction", bool),
        "NCORES": ("n_cores", int),
        "FIND_FILE": ("find_file", bool),
        "INIT_JWST": ("execute_jwst_initialisation", bool),
        "GEN_PSF": ("generate_psf", bool),
        "GEN_RUN": ("generate_run", bool),
        "GEN_REGION": ("generate_region", bool),
        "AST_RECOVER": ("ast_recover", bool),
        "AST_AUTOSAVE": ("ast_auto_save", int),
        "AST_NO_BGD": ("ast_no_background", bool),
        "AST_NO_PSF": ("ast_no_psf_phot", bool),
        "MATCH_BAND": ("do_band_processing", bool),
        "MATCH_CASCADE": ("do_cascade", bool),
        "MATCH_DITHER": ("use_dither", bool),
        "MATCH_EXACT": ("exact_match", bool),
        "MATCH_FULL": ("full_run", bool),
        "MATCH_GENERIC": ("generic_mode", bool),
        "MATCH_ERR_COL": ("error_col", str),
        "MATCH_MASK_EVAL": ("mask_eval", str),
        "PLOT_TEST": ("test_mode", bool),
        "PLOT_DARK": ("dark_mode", bool),
        "PLOT_INSPECT": ("inspect_parameter", str),
        "PLOT_STYLE": ("plot_style", str),
        "PARAM": ("param_tag", str),
    }

    def __init__(self) -> None:
        self._frozen: bool = False

        # high level stuff
        self._show_help: bool = False
        self._show_ast_help: bool = False
        self._verbose_logs: bool = False
        self._show_version: bool = False

        # main actions
        self._do_aperture_photometry: bool = False
        self._do_bgd_estimate: bool = False
        self._do_star_detection: bool = False
        self._do_source_geometry: bool = False
        self._do_matching: bool = False
        self._do_photometry_routine: bool = False
        self._do_bgd_subtraction: bool = False

        # other actions
        self._generate_psf: bool = False
        self._generate_run: bool = False
        self._generate_region: bool = False
        self._generate_local_param_file: bool = False
        self._execute_jwst_initialisation = False

        # updates an old version param file.
        self._update_param: bool = False

        # file parameters
        self._param_file: str | None = None
        self._ap_file: str | None = None
        self._background_file: str | None = None
        self._find_file: bool = True
        self._region_file: str | None = None

        # multiprocessing params (assumes to use 1 core to begin with)
        self._n_cores: int = 1

        # artificial stars params
        self._do_artificial_star_test: bool = False
        self._do_artificial_star_test_results: bool = False
        self._ast_out_tables: list[Table] | None = None
        self._ast_recover: bool = False
        self._ast_auto_save: int = 100
        self._ast_no_background: bool = True
        self._ast_no_psf_phot: bool = True
        self._save_added_image: bool = False
        self._save_added_image_path: str = ""
        self._ast_seed: int | None = None

        # states saved in config as only access route into the runs for
        # optimization purposes
        self._ast_loader: np.ndarray | None = None
        self._ast_test_index: int = 0
        self._ast_psf: np.ndarray | None = None
        self._ast_load_psf: bool = True
        self._ast_add_stars: bool = False

        # matching params
        self._do_band_processing: bool = False
        self._do_cascade: bool = False
        self._use_dither: bool = False
        self._exact_match: bool = False
        self._full_run: bool = False
        self._generic_mode: bool = False
        self._show_match_help: bool = False
        self._band_deprecated: bool = False
        self._error_col: str = TableColumn.E_FLUX
        self._mask_eval: str | None = None

        # plot params
        self._show_plot_help: bool = False
        self._test_mode: bool = False
        self._dark_frame_correction = False
        self._inspect_parameter: str | None = None
        self._plot_style: str | None = None

        # param file defaults. These constants do not have justifications yet.
        self._output_file: str | None = os.getenv(STARBUG_DATA_DIR)
        self._hdu_name: str = SCI
        self._filter: str | None = None
        self._full_width_half_max: float = -1.0
        self._sigma_sky: float = 2.0
        self._sigma_source: float = 5.0
        self._do_bgd_2d: bool = True
        self._do_convolution: bool = True
        self._clean_sources: bool = True
        self._sharp_cutoff_low: float = 0.4
        self._sharp_cutoff_high: float = 0.9
        self._round1_cutoff_high: float = 1.0
        self._round2_cutoff_high: float = 1.0
        self._smooth_low: float = 0.0
        self._smooth_high: float = 1.0
        self._ricker_wavelet_radius: float = 1.0
        self._aperture_phot_radius: float = 1.5
        self._encircled_energy_fraction: float = -1.0
        self._sky_annulus_inner_radius: float = 3.0
        self._sky_annulus_outer_radius: float = 4.5
        self._ap_corr_file_override: str | None = None
        self._bgd_radius: float = 0.0
        self._profile_scaling_factor: float = 1.0
        self._profile_slope: float = 0.5
        self._background_box_size: int = 2
        self._bgd_check_file: str | None = None
        self._psf_file_override: str = DEFAULT_PSF_FILE_NAME
        self._use_wcs_values: bool = True
        self._zero_point_magnitude: float = 8.9
        self._critical_separation: float = 1.0
        self._force_centroid_position: bool = False
        self._centroid_delta_threshold: float = -1.0
        self._max_xy_deviation: str = '3.0'
        self._psf_fit_size: int = -1
        self._generate_residual_image: bool = False
        self._calculate_crowding_metric: bool = True
        self._match_threshold_arc_sec: str = "0.1"
        self._extra_match_columns: str | None = None
        self._exposure_count_threshold: int = -1
        self._signal_to_noise_threshold: float = -1.0
        self._bridge_band_column: str | None = None
        self._artificial_star_tests_count: int = 100
        self._stars_per_artificial_test: int = 10
        self._sub_image_crop_size: int = 500
        self._test_magnitude_bright_limit: int = DEFAULT_MAX_MAG
        self._test_magnitude_faint_limit: int = DEFAULT_MIN_MAG
        self._ast_plot_filename: str | None = None
        self._region_colour: str = DEFAULT_COLOUR
        self._region_scale: bool = True
        self._region_radius: int = 3
        self._region_x_column_name: str = TableColumn.RA
        self._region_y_column_name: str = TableColumn.DEC
        self._region_uses_wcs: bool = True
        self._param_tag: str = STAR_BUG_PARAMS

        # generate psf variables
        self._detector_name: str | None = None

        # target images
        self._fits_images: list[str] = []

    @classmethod
    def _generate_get_opt_definitions(cls, param_map) -> Tuple[str, list[str]]:
        # noinspection SpellCheckingInspection
        """
        Natively inspects the configuration structure to construct perfectly
        formatted short-option strings and long-option arrays for getopt
        dynamically.

        :return: the inputs to gnu_getopt.
        """
        short_opts_compiled = ""
        long_opts_compiled = []

        for (short_flag, long_flag, data_type) in param_map.keys():
            suffix = "" if data_type is bool else "="
            long_opts_compiled.append(f"{long_flag}{suffix}")

            if short_flag:
                short_suffix = "" if data_type is bool else ":"
                short_opts_compiled += f"{short_flag}{short_suffix}"

        return short_opts_compiled, long_opts_compiled

    @classmethod
    def generate_main_get_opt_definitions(cls) -> Tuple[str, list[str]]:
        # noinspection SpellCheckingInspection
        """
        Natively inspects the configuration structure to construct perfectly
        formatted short-option strings and long-option arrays for getopt
        dynamically for main.

        :return: the inputs to gnu_getopt.
        """
        return cls._generate_get_opt_definitions(cls.MAIN_FLAG_MAP)

    @classmethod
    def generate_ast_get_opt_definitions(cls) -> Tuple[str, list[str]]:
        # noinspection SpellCheckingInspection
        """
        Natively inspects the configuration structure to construct perfectly
        formatted short-option strings and long-option arrays for getopt
        dynamically for artificial stars.

        :return: the inputs to gnu_getopt.
        """
        return cls._generate_get_opt_definitions(cls.AST_FLAG_MAP)

    @classmethod
    def generate_match_get_opt_definitions(cls) -> Tuple[str, list[str]]:
        # noinspection SpellCheckingInspection
        """
        Natively inspects the configuration structure to construct perfectly
        formatted short-option strings and long-option arrays for getopt
        dynamically for match.

        :return: the inputs to gnu_getopt.
        """
        return cls._generate_get_opt_definitions(cls.MATCH_FLAG_MAP)

    @classmethod
    def generate_plot_get_opt_definitions(cls) -> Tuple[str, list[str]]:
        # noinspection SpellCheckingInspection
        """
        Natively inspects the configuration structure to construct perfectly
        formatted short-option strings and long-option arrays for getopt
        dynamically for plot.

        :return: the inputs to gnu_getopt.
        """
        return cls._generate_get_opt_definitions(cls.PLOT_FLAG_MAP)

    @staticmethod
    def parse_value(val_str: str) -> int | float | str:
        # Clean up tracking whitespaces while they are guaranteed to be
        # strings
        raw_value = val_str.strip()

        # Default fallback type is the cleaned string itself
        value: int | float | str = raw_value

        # Attempt numeric type conversions safely
        try:
            if '.' in raw_value:
                value = float(raw_value)
            else:
                value = int(raw_value)
        except (ValueError, AttributeError, TypeError):
            # If conversion fails, value remains a string
            pass
        return value

    @staticmethod
    def parse_param(line: str) -> Dict[str, int | float | str]:
        """
        Parse a parameter line
        :param line: the line to parse
        :type line: str
        :return: the parsed params from the line
        :rtype: Dict[str, int | float | str]
        """
        param: Dict[str, int | float | str] = {}

        # Guard against empty lines or comment lines
        if line and line[0] not in StarBugMainConfig.EXCLUDES:
            key_str: str = ""
            val_str: str = ""

            # Explicitly parse out raw string substrings first
            if "//" in line and line[0] != "/":
                parsed = parse("{}={}//{}", line)
                if parsed:
                    key_str, val_str, _ = parsed
            else:
                parsed = parse("{}={}", line)
                if parsed:
                    key_str, val_str = parsed

            # Clean up tracking whitespaces while they are guaranteed to be
            # strings
            key = key_str.strip()

            value: int | float | str
            value = StarBugMainConfig.parse_value(val_str)

            # Special case environmental variables expansions for paths
            if (key in (HeaderTags.OUTPUT, AP_FILE, BGD_FILE, PSF_FILE)
                    and isinstance(value, str)):
                value = os.path.expandvars(value)

            if value == "False":
                raise Exception("")
            param[key] = value

        return param

    @staticmethod
    def load_params(f_name) -> 'StarBugMainConfig':
        """
        Convert a parameter file into a dictionary of options

        :param f_name: path/to/file.param
        :type f_name: str or None
        :return: dictionary of options
        :rtype: dict of string, string
        """
        config: StarBugMainConfig = StarBugMainConfig()
        if f_name is None:
            return config

        if os.path.exists(f_name):
            with open(f_name, "r") as fp:
                for line in fp.readlines():
                    config.update(StarBugMainConfig.parse_param(line))
        else:
            p_error("config file \"%s\" does not exist\n. Using "
                    "default config instead" % f_name)
        return config

    def populate_params(
            self, argv: list[str], short_definition: str,
            long_definition: list[str],
            param_map: dict[tuple[str | None, str, Any], str]) -> None:
        """
        populates the config from command line requests
        :param argv: the command line
        :param short_definition: the short lists of command lines
        :param long_definition: the large list of command lines
        :param param_map: mapping between command line and property.
        :return: None
        """
        opts: list[tuple[str, str]]
        args: list[str]
        opts, args = getopt.gnu_getopt(
            argv, short_definition, long_definition)

        for opt, opt_arg in opts:
            # strip down to raw option label text
            clean_opt = opt.lstrip('-')

            for (short_flag, long_flag, data_type), property_name in (
                    param_map.items()):
                if clean_opt in (short_flag, long_flag):
                    if data_type is bool:
                        setattr(self, property_name, True)
                    else:
                        # Casts string inputs to explicit types for
                        # variables
                        setattr(self, property_name, data_type(opt_arg))
                    break

        # add the training args as target image files
        self.fits_images = args

    def got_valid_psf_generation_params(self) -> bool:
        """
        returns if the config has parameters set correctly to execute psf
        generation
        :return: bool if the params are set away from invalid defaults.
        :rtype: bool
        """
        return (self._filter != "" and self._detector_name != ""
                and self._psf_fit_size != -1)

    def use_main_one_time_runs(self) -> bool:
        """
        check for any of the one-off runs.
        :return: bool if there is one time runs to run.
        :rtype: bool
        """
        return (
            self._show_help or self._update_param or
            self._execute_jwst_initialisation or self._generate_psf or
            self._generate_run or self._generate_region or
            self._generate_local_param_file or self._show_version)

    def use_ast_one_time_runs(self) -> bool:
        """
        check for any of the one-off runs.
        :return: bool if there is one time runs to run.
        :rtype: bool
        """
        return self._show_ast_help or self._ast_recover

    def generate_default_param_file_text(self, version_str: str) -> str:
        """
        Dynamically constructs the entire default configuration string template
        using the active internal variable defaults directly via the named
        constant template.
        """
        format_dictionary: dict[str, str] = {}
        for key, (prop, target_type) in self.MAIN_PARAM_FILE_MAP.items():
            val = getattr(self, prop)
            if target_type is bool:
                format_dictionary[key] = "1" if val else "0"
            else:
                format_dictionary[key] = "" if val is None else str(val)

        # populate default param file.
        default_template: Path = (
            Path(__file__).parent / "default_param_file.txt")
        template_str = default_template.read_text(encoding="utf-8")
        return template_str.format(
            version_str=version_str, **format_dictionary)

    def do_generate_local_param_file(self) -> None:
        """
        writes a local param file based off the state of this config.
        :return: None
        """
        # handle output file
        path: str = "starbug.param"
        if self._output_file is not None:
            path = os.path.join(self._output_file, "starbug.param")

        # generate new param file
        with open(path, "w") as fp:
            fp.write(self.generate_default_param_file_text(get_version()))

    def update(self, update_values: dict[str, str | int | float]) -> None:
        """
        updates the config with new values extracted from a param file.
        :param update_values: the updated values
        :return: None
        """
        for key, raw_value in update_values.items():
            if key not in self.MAIN_PARAM_FILE_MAP.keys():
                raise TypeError(
                    f"Param {key} no longer works within Starbug2. Please "
                    f"execute starbug2 --update-param")
            property_name, target_type = self.MAIN_PARAM_FILE_MAP[key]

            # process raw value
            if raw_value == "" or raw_value is None:
                setattr(self, property_name, None)
                continue

            if target_type is bool:
                # Converts parameter flags like 0 or 1 integers to standard
                # Booleans
                setattr(self, property_name, bool(int(raw_value)))
            else:
                setattr(self, property_name, target_type(raw_value))

    def _normalize_threshold(
            self, threshold: float | int | np.ndarray | list | Quantity) -> (
            None | np.ndarray | Quantity):
        # noinspection SpellCheckingInspection
        """
        Normalises threshold inputs to ensure they possess the 'arcsec' unit.
        - Unitless Quantities -> scaled to arcsec
        - Floats/Ints -> converted to arcsec Quantity
        - Lists/Arrays of floats -> converted to an object array of arcsec
                                    Quantities
        """
        if threshold is None:
            return None

        # Handle standard Lists or NumPy arrays
        if isinstance(threshold, (list, np.ndarray)):
            # Recursively normalise each element, keeping them as individual
            # object array slots
            return np.array(
                [self._normalize_threshold(t) for t in threshold],
                dtype=object)

        # Handle Astropy Quantity instances
        if isinstance(threshold, Quantity):
            # Check if the Quantity has no physical units (dimensionless)
            if threshold.unit == units.dimensionless_unscaled:
                return threshold.value * units.arcsec
            return threshold

        # Handle raw primitive Python numeric scalars (ints, floats)
        if isinstance(threshold, (int, float, np.number)):
            return threshold * units.arcsec
        return None

    def freeze(self) -> None:
        """
        locks the class from being editable
        :return: None
        """
        self._frozen = True

    def unfreeze(self) -> None:
        """
        unlocks the class to become editable.
        NOTE: if frozen previously. this really shouldn't happen.
        :return: None
        """
        self._frozen = False

    # ==========================================
    # these 2 methods are here to ensure the access to threshold is in
    # a quantity and an array of quantity as needed
    # ==========================================
    @property
    def match_threshold_arc_sec_as_an_arc_sec(self) -> units.Quantity:
        """
        builds threshold as an arc second Quantity.
        :return: the threshold as a quantity
        :rtype: units.Quantity
        """
        threshold: float = float(self._match_threshold_arc_sec)
        normalised_threshold: None | np.ndarray | Quantity = (
            self._normalize_threshold(threshold))
        assert isinstance(normalised_threshold, Quantity)
        return normalised_threshold

    @property
    def match_threshold_arc_sec_as_an_array(self) -> np.ndarray:
        """
        builds threshold as an array of quantity
        :return: the threshold as an array of quantity
        :rtype: ndarray[Quantity]
        """
        threshold_array: np.ndarray = np.array(
            self._match_threshold_arc_sec.split(','), float)
        normalised_threshold: None | np.ndarray | Quantity = (
            self._normalize_threshold(threshold_array))
        assert isinstance(normalised_threshold, np.ndarray)
        return normalised_threshold

    # ==========================================
    # BELOW HERE ARE GETTERS AND SETTERS FOR EVERYTHING
    # ==========================================

    # ==========================================
    # HIGH LEVEL STUFF
    # ==========================================

    @property
    def show_help(self) -> bool:
        return self._show_help

    @show_help.setter
    def show_help(self, value: bool) -> None:
        self._show_help = value

    @property
    def verbose_logs(self) -> bool:
        return self._verbose_logs

    @verbose_logs.setter
    def verbose_logs(self, value: bool) -> None:
        self._verbose_logs = value

    @property
    def show_version(self) -> bool:
        return self._show_version

    @show_version.setter
    def show_version(self, value: bool) -> None:
        self._show_version = value

    # ==========================================
    # MAIN ACTIONS
    # ==========================================

    @property
    def do_aperture_photometry(self) -> bool:
        return self._do_aperture_photometry

    @do_aperture_photometry.setter
    def do_aperture_photometry(self, value: bool) -> None:
        self._do_aperture_photometry = value

    @property
    def do_bgd_estimate(self) -> bool:
        return self._do_bgd_estimate

    @do_bgd_estimate.setter
    def do_bgd_estimate(self, value: bool) -> None:
        self._do_bgd_estimate = value

    @property
    def do_star_detection(self) -> bool:
        return self._do_star_detection

    @do_star_detection.setter
    def do_star_detection(self, value: bool) -> None:
        self._do_star_detection = value
        self._do_aperture_photometry = value

    @property
    def do_source_geometry(self) -> bool:
        return self._do_source_geometry

    @do_source_geometry.setter
    def do_source_geometry(self, value: bool) -> None:
        self._do_source_geometry = value

    @property
    def do_matching(self) -> bool:
        return self._do_matching

    @do_matching.setter
    def do_matching(self, value: bool) -> None:
        self._do_matching = value

    @property
    def do_photometry_routine(self) -> bool:
        return self._do_photometry_routine

    @do_photometry_routine.setter
    def do_photometry_routine(self, value: bool) -> None:
        self._do_photometry_routine = value

    @property
    def do_bgd_subtraction(self) -> bool:
        return self._do_bgd_subtraction

    @do_bgd_subtraction.setter
    def do_bgd_subtraction(self, value: bool) -> None:
        self._do_bgd_subtraction = value

    # ==========================================
    # OTHER ACTIONS
    # ==========================================

    @property
    def generate_psf(self) -> bool:
        return self._generate_psf

    @generate_psf.setter
    def generate_psf(self, value: bool) -> None:
        self._generate_psf = value

    @property
    def generate_run(self) -> bool:
        return self._generate_run

    @generate_run.setter
    def generate_run(self, value: bool) -> None:
        self._generate_run = value

    @property
    def generate_region(self) -> bool:
        return self._generate_region

    @generate_region.setter
    def generate_region(self, value: bool) -> None:
        self._generate_region = value

    @property
    def generate_local_param_file(self) -> bool:
        return self._generate_local_param_file

    @generate_local_param_file.setter
    def generate_local_param_file(self, value: bool) -> None:
        self._generate_local_param_file = value

    @property
    def update_param(self) -> bool:
        return self._update_param

    @update_param.setter
    def update_param(self, value: bool) -> None:
        self._update_param = value

    # ==========================================
    # PARAMETERS
    # ==========================================

    @property
    def param_file(self) -> str | None:
        return self._param_file

    @param_file.setter
    def param_file(self, value: str | None) -> None:
        self._param_file = value

    @property
    def ap_file(self) -> str | None:
        return self._ap_file

    @ap_file.setter
    def ap_file(self, value: str | None) -> None:
        if value is None or os.path.exists(value):
            self._ap_file = value
        else:
            p_error("AP_FILE \"%s\" does not exist\n" % value)

    @property
    def background_file(self) -> str | None:
        return self._background_file

    @background_file.setter
    def background_file(self, value: str | None) -> None:
        if value is None or os.path.exists(value):
            self._background_file = value
        else:
            p_error("BGD_FILE \"%s\" does not exist\n" % value)

    @property
    def find_file(self) -> bool:
        return self._find_file

    @find_file.setter
    def find_file(self, value: bool) -> None:
        self._find_file = value

    @property
    def n_cores(self) -> int:
        return self._n_cores

    @n_cores.setter
    def n_cores(self, value: int) -> None:
        self._n_cores = value

    # ==========================================
    # DYNAMIC PARAM FILE PROPERTIES (MONSTER EXTRAS)
    # ==========================================

    @property
    def output_file(self) -> str | None:
        return self._output_file

    @output_file.setter
    def output_file(self, value: str) -> None:
        self._output_file = value

    @property
    def hdu_name(self) -> str | int | float:
        return StarBugMainConfig.parse_value(self._hdu_name)

    @hdu_name.setter
    def hdu_name(self, value: str) -> None:
        self._hdu_name = value

    @property
    def custom_filter(self) -> str | None:
        return self._filter

    @custom_filter.setter
    def custom_filter(self, value: str) -> None:
        # added warning if we're planning on using F150W2 filter, as currently
        # issue arises that we've had to 1/2 the resolution to allow it to
        # pass init. see https://github.com/alan-stokes/starbug2/issues/2
        # for more details.
        self._filter = value

        if self._filter == PROBLEMATIC_FILTER_ID:
            warn(PROBLEMATIC_FILTER_WARNING)

    @property
    def full_width_half_max(self) -> float:
        return self._full_width_half_max

    @full_width_half_max.setter
    def full_width_half_max(self, value: float) -> None:
        self._full_width_half_max = value

    def full_width_half_max_with_filter(
            self, filter_struct: FilterStruct | None) -> float:
        """
        figure the full_width_half_max from the filter struct
        :param filter_struct: the filter struct.
        :type filter_struct: FilterStruct
        :return: the full_width_half_max
        :rtype: float
        """
        if self._full_width_half_max > 0:
            return self._full_width_half_max
        elif filter_struct:
            return filter_struct.full_width_half_max
        else:
            return DEFAULT_FULL_WIDTH_HALF_MAX

    @property
    def sigma_sky(self) -> float:
        return self._sigma_sky

    @sigma_sky.setter
    def sigma_sky(self, value: float) -> None:
        self._sigma_sky = value

    @property
    def sigma_source(self) -> float:
        return self._sigma_source

    @sigma_source.setter
    def sigma_source(self, value: float) -> None:
        self._sigma_source = value

    @property
    def do_bgd_2d(self) -> bool:
        return self._do_bgd_2d

    @do_bgd_2d.setter
    def do_bgd_2d(self, value: bool) -> None:
        self._do_bgd_2d = value

    @property
    def do_convolution(self) -> bool:
        return self._do_convolution

    @do_convolution.setter
    def do_convolution(self, value: bool) -> None:
        self._do_convolution = value

    @property
    def clean_sources(self) -> bool:
        return self._clean_sources

    @clean_sources.setter
    def clean_sources(self, value: bool) -> None:
        self._clean_sources = value

    @property
    def sharp_cutoff_low(self) -> float:
        return self._sharp_cutoff_low

    @sharp_cutoff_low.setter
    def sharp_cutoff_low(self, value: float) -> None:
        self._sharp_cutoff_low = value

    @property
    def sharp_cutoff_high(self) -> float:
        return self._sharp_cutoff_high

    @sharp_cutoff_high.setter
    def sharp_cutoff_high(self, value: float) -> None:
        self._sharp_cutoff_high = value

    @property
    def round1_cutoff_high(self) -> float:
        return self._round1_cutoff_high

    @round1_cutoff_high.setter
    def round1_cutoff_high(self, value: float) -> None:
        self._round1_cutoff_high = value

    @property
    def round2_cutoff_high(self) -> float:
        return self._round2_cutoff_high

    @round2_cutoff_high.setter
    def round2_cutoff_high(self, value: float) -> None:
        self._round2_cutoff_high = value

    @property
    def smooth_low(self) -> float:
        return self._smooth_low

    @smooth_low.setter
    def smooth_low(self, value: float) -> None:
        self._smooth_low = value

    @property
    def smooth_high(self) -> float:
        return self._smooth_high

    @smooth_high.setter
    def smooth_high(self, value: float) -> None:
        self._smooth_high = value

    @property
    def ricker_wavelet_radius(self) -> float:
        return self._ricker_wavelet_radius

    @ricker_wavelet_radius.setter
    def ricker_wavelet_radius(self, value: float) -> None:
        self._ricker_wavelet_radius = value

    @property
    def aperture_phot_radius(self) -> float:
        return self._aperture_phot_radius

    @aperture_phot_radius.setter
    def aperture_phot_radius(self, value: float) -> None:
        self._aperture_phot_radius = value

    @property
    def encircled_energy_fraction(self) -> float:
        return self._encircled_energy_fraction

    @encircled_energy_fraction.setter
    def encircled_energy_fraction(self, value: float) -> None:
        self._encircled_energy_fraction = value

    @property
    def sky_annulus_inner_radius(self) -> float:
        return self._sky_annulus_inner_radius

    @sky_annulus_inner_radius.setter
    def sky_annulus_inner_radius(self, value: float) -> None:
        self._sky_annulus_inner_radius = value

    @property
    def sky_annulus_outer_radius(self) -> float:
        return self._sky_annulus_outer_radius

    @sky_annulus_outer_radius.setter
    def sky_annulus_outer_radius(self, value: float) -> None:
        self._sky_annulus_outer_radius = value

    @property
    def ap_corr_file_override(self) -> str | None:
        return self._ap_corr_file_override

    @ap_corr_file_override.setter
    def ap_corr_file_override(self, value: str) -> None:
        self._ap_corr_file_override = value

    @property
    def bgd_radius(self) -> float:
        return self._bgd_radius

    @bgd_radius.setter
    def bgd_radius(self, value: float) -> None:
        self._bgd_radius = value

    @property
    def profile_scaling_factor(self) -> float:
        return self._profile_scaling_factor

    @profile_scaling_factor.setter
    def profile_scaling_factor(self, value: float) -> None:
        self._profile_scaling_factor = value

    @property
    def profile_slope(self) -> float:
        return self._profile_slope

    @profile_slope.setter
    def profile_slope(self, value: float) -> None:
        self._profile_slope = value

    @property
    def background_box_size(self) -> int:
        return self._background_box_size

    @background_box_size.setter
    def background_box_size(self, value: int) -> None:
        self._background_box_size = value

    @property
    def bgd_check_file(self) -> str | None:
        return self._bgd_check_file

    @bgd_check_file.setter
    def bgd_check_file(self, value: str) -> None:
        self._bgd_check_file = value

    @property
    def psf_file_override(self) -> str:
        return self._psf_file_override

    @psf_file_override.setter
    def psf_file_override(self, value: str) -> None:
        self._psf_file_override = value

    @property
    def use_wcs_values(self) -> bool:
        return self._use_wcs_values

    @use_wcs_values.setter
    def use_wcs_values(self, value: bool) -> None:
        self._use_wcs_values = value

    @property
    def zero_point_magnitude(self) -> float:
        return self._zero_point_magnitude

    @zero_point_magnitude.setter
    def zero_point_magnitude(self, value: float) -> None:
        self._zero_point_magnitude = value

    @property
    def critical_separation(self) -> float:
        return self._critical_separation

    @critical_separation.setter
    def critical_separation(self, value: float) -> None:
        self._critical_separation = value

    @property
    def force_centroid_position(self) -> bool:
        return self._force_centroid_position

    @force_centroid_position.setter
    def force_centroid_position(self, value: bool) -> None:
        self._force_centroid_position = value

    @property
    def max_xy_deviation(self) -> str:
        return self._max_xy_deviation

    @max_xy_deviation.setter
    def max_xy_deviation(self, value: str) -> None:
        self._max_xy_deviation = value

    @property
    def psf_fit_size(self) -> int:
        return self._psf_fit_size

    @psf_fit_size.setter
    def psf_fit_size(self, value: int) -> None:
        self._psf_fit_size = value

    @property
    def generate_residual_image(self) -> bool:
        return self._generate_residual_image

    @generate_residual_image.setter
    def generate_residual_image(self, value: bool) -> None:
        self._generate_residual_image = value

    @property
    def calculate_crowding_metric(self) -> bool:
        return self._calculate_crowding_metric

    @calculate_crowding_metric.setter
    def calculate_crowding_metric(self, value: bool) -> None:
        self._calculate_crowding_metric = value

    @property
    def match_threshold_arc_sec(self) -> str:
        return self._match_threshold_arc_sec

    @match_threshold_arc_sec.setter
    def match_threshold_arc_sec(self, value: str) -> None:
        self._match_threshold_arc_sec = value

    @property
    def extra_match_columns(self) -> str:
        if self._extra_match_columns is None:
            return ""
        return self._extra_match_columns

    @extra_match_columns.setter
    def extra_match_columns(self, value: str) -> None:
        self._extra_match_columns = value

    @property
    def exposure_count_threshold(self) -> int:
        return self._exposure_count_threshold

    @exposure_count_threshold.setter
    def exposure_count_threshold(self, value: int) -> None:
        self._exposure_count_threshold = value

    @property
    def bridge_band_column(self) -> str:
        if self._bridge_band_column is None:
            return ""
        return self._bridge_band_column

    @bridge_band_column.setter
    def bridge_band_column(self, value: str) -> None:
        self._bridge_band_column = value

    @property
    def artificial_star_tests_count(self) -> int:
        return self._artificial_star_tests_count

    @artificial_star_tests_count.setter
    def artificial_star_tests_count(self, value: int) -> None:
        self._artificial_star_tests_count = value

    @property
    def stars_per_artificial_test(self) -> int:
        return self._stars_per_artificial_test

    @stars_per_artificial_test.setter
    def stars_per_artificial_test(self, value: int) -> None:
        self._stars_per_artificial_test = value

    @property
    def sub_image_crop_size(self) -> int:
        return self._sub_image_crop_size

    @sub_image_crop_size.setter
    def sub_image_crop_size(self, value: int) -> None:
        self._sub_image_crop_size = value

    @property
    def test_magnitude_bright_limit(self) -> int:
        return self._test_magnitude_bright_limit

    @test_magnitude_bright_limit.setter
    def test_magnitude_bright_limit(self, value: int) -> None:
        self._test_magnitude_bright_limit = value

    @property
    def test_magnitude_faint_limit(self) -> int:
        return self._test_magnitude_faint_limit

    @test_magnitude_faint_limit.setter
    def test_magnitude_faint_limit(self, value: int) -> None:
        self._test_magnitude_faint_limit = value

    @property
    def ast_plot_filename(self) -> str | None:
        return self._ast_plot_filename

    @ast_plot_filename.setter
    def ast_plot_filename(self, value: str) -> None:
        self._ast_plot_filename = value

    @property
    def region_colour(self) -> str:
        return self._region_colour

    @region_colour.setter
    def region_colour(self, value: str) -> None:
        self._region_colour = value

    @property
    def region_scale(self) -> bool:
        return self._region_scale

    @region_scale.setter
    def region_scale(self, value: bool) -> None:
        self._region_scale = value

    @property
    def region_radius(self) -> int:
        return self._region_radius

    @region_radius.setter
    def region_radius(self, value: int) -> None:
        self._region_radius = value

    @property
    def region_x_column_name(self) -> str:
        return self._region_x_column_name

    @region_x_column_name.setter
    def region_x_column_name(self, value: str) -> None:
        self._region_x_column_name = value

    @property
    def region_y_column_name(self) -> str:
        return self._region_y_column_name

    @region_y_column_name.setter
    def region_y_column_name(self, value: str) -> None:
        self._region_y_column_name = value

    @property
    def region_uses_wcs(self) -> bool:
        return self._region_uses_wcs

    @region_uses_wcs.setter
    def region_uses_wcs(self, value: bool) -> None:
        self._region_uses_wcs = value

    @property
    def execute_jwst_initialisation(self) -> bool:
        return self._execute_jwst_initialisation

    @execute_jwst_initialisation.setter
    def execute_jwst_initialisation(self, value: bool) -> None:
        self._execute_jwst_initialisation = value

    @property
    def detector_name(self) -> str | None:
        return self._detector_name

    @detector_name.setter
    def detector_name(self, value: str) -> None:
        self._detector_name = value

    @property
    def region_file(self) -> str | None:
        return self._region_file

    @region_file.setter
    def region_file(self, value: str) -> None:
        self._region_file = value

    @property
    def fits_images(self) -> list[str]:
        return self._fits_images

    @fits_images.setter
    def fits_images(self, values: list[str]) -> None:
        self._fits_images = values

    # when in match, it is no longer image data. but table data. so utilise
    # this method for clarity of user readability.
    @property
    def fits_table(self) -> list[str]:
        return self._fits_images

    @fits_table.setter
    def fits_table(self, values: list[str]) -> None:
        self._fits_images = values

    @property
    def param_tag(self) -> str:
        return self._param_tag

    @param_tag.setter
    def param_tag(self, value) -> None:
        self._param_tag = value

    # ===============================
    # AST properties
    # ===============================

    @property
    def ast_add_stars(self) -> bool:
        return self._ast_add_stars

    @ast_add_stars.setter
    def ast_add_stars(self, value) -> None:
        self._ast_add_stars = value

    @property
    def ast_load_psf(self) -> bool:
        return self._ast_load_psf

    @ast_load_psf.setter
    def ast_load_psf(self, value: bool) -> None:
        self._ast_load_psf = value

    @property
    def ast_psf(self) -> np.ndarray | None:
        return self._ast_psf

    @ast_psf.setter
    def ast_psf(self, value: np.ndarray | None) -> None:
        self._ast_psf = value

    @property
    def do_artificial_star_test(self) -> bool:
        return self._do_artificial_star_test

    @do_artificial_star_test.setter
    def do_artificial_star_test(self, value: bool) -> None:
        self._do_artificial_star_test = value

    @property
    def ast_loader(self) -> np.ndarray | None:
        return self._ast_loader

    @ast_loader.setter
    def ast_loader(self, value: np.ndarray) -> None:
        self._ast_loader = value

    @property
    def ast_test_index(self) -> int:
        return self._ast_test_index

    @ast_test_index.setter
    def ast_test_index(self, value: int) -> None:
        self._ast_test_index = value

    @property
    def ast_seed(self) -> int | None:
        return self._ast_seed

    @ast_seed.setter
    def ast_seed(self, value: int | None) -> None:
        self._ast_seed = value

    @property
    def ast_save_added_image_path(self) -> str:
        return self._save_added_image_path

    @ast_save_added_image_path.setter
    def ast_save_added_image_path(self, value: str) -> None:
        self._save_added_image_path = value

    @property
    def ast_save_added_image(self) -> bool:
        return self._save_added_image

    @ast_save_added_image.setter
    def ast_save_added_image(self, value: bool) -> None:
        self._save_added_image = value

    @property
    def show_ast_help(self) -> bool:
        return self._show_ast_help

    @show_ast_help.setter
    def show_ast_help(self, value: bool) -> None:
        self._show_ast_help = value

    @property
    def ast_recover(self) -> bool:
        return self._ast_recover

    @ast_recover.setter
    def ast_recover(self, value: bool) -> None:
        self._ast_recover = value

    @property
    def ast_auto_save(self) -> int:
        return self._ast_auto_save

    @ast_auto_save.setter
    def ast_auto_save(self, value: int) -> None:
        self._ast_auto_save = value

    @property
    def ast_no_background(self) -> bool:
        return self._ast_no_background

    @ast_no_background.setter
    def ast_no_background(self, value: bool) -> None:
        self._ast_no_background = value

    @property
    def ast_no_psf_phot(self) -> bool:
        return self._ast_no_psf_phot

    @ast_no_psf_phot.setter
    def ast_no_psf_phot(self, value: bool) -> None:
        self._ast_no_psf_phot = value

    @property
    def do_artificial_star_test_results(self) -> bool:
        return self._do_artificial_star_test_results

    @do_artificial_star_test_results.setter
    def do_artificial_star_test_results(self, value: bool) -> None:
        self._do_artificial_star_test_results = value

    @property
    def ast_out_tables(self) -> list[Table] | None:
        return self._ast_out_tables

    @ast_out_tables.setter
    def ast_out_tables(self, value: list[Table]) -> None:
        self._ast_out_tables = value

    # =======================================
    # matching properties
    # =======================================

    @property
    def do_band_processing(self) -> bool:
        return self._do_band_processing

    @do_band_processing.setter
    def do_band_processing(self, value: bool) -> None:
        self._do_band_processing = value

    @property
    def do_cascade(self) -> bool:
        return self._do_cascade

    @do_cascade.setter
    def do_cascade(self, value: bool) -> None:
        self._do_cascade = value

    @property
    def use_dither(self) -> bool:
        return self._use_dither

    @use_dither.setter
    def use_dither(self, value: bool) -> None:
        self._use_dither = value

    @property
    def exact_match(self) -> bool:
        return self._exact_match

    @exact_match.setter
    def exact_match(self, value: bool) -> None:
        self._exact_match = value

    @property
    def full_run(self) -> bool:
        return self._full_run

    @full_run.setter
    def full_run(self, value: bool) -> None:
        self._full_run = value

    @property
    def generic_mode(self) -> bool:
        return self._generic_mode

    @generic_mode.setter
    def generic_mode(self, value: bool) -> None:
        self._generic_mode = value

    @property
    def show_match_help(self) -> bool:
        return self._show_match_help

    @show_match_help.setter
    def show_match_help(self, value: bool) -> None:
        self._show_match_help = value

    @property
    def band_deprecated(self) -> bool:
        return self._band_deprecated

    @band_deprecated.setter
    def band_deprecated(self, value: bool) -> None:
        self._band_deprecated = value

    @property
    def error_col(self) -> str:
        return self._error_col

    @error_col.setter
    def error_col(self, value: str) -> None:
        self._error_col = value

    @property
    def mask_eval(self) -> str | None:
        return self._mask_eval

    @mask_eval.setter
    def mask_eval(self, value: str | None) -> None:
        self._mask_eval = value

    # ==============================
    # plot getters and setters
    # ==============================

    @property
    def show_plot_help(self) -> bool:
        return self._show_plot_help

    @show_plot_help.setter
    def show_plot_help(self, value: bool) -> None:
        self._show_plot_help = value

    @property
    def test_mode(self) -> bool:
        return self._test_mode

    @test_mode.setter
    def test_mode(self, value: bool) -> None:
        self._test_mode = value

    @property
    def dark_mode(self) -> bool:
        return self._dark_frame_correction

    @dark_mode.setter
    def dark_mode(self, value: bool) -> None:
        self._dark_frame_correction = value

    @property
    def inspect_parameter(self) -> str | None:
        return self._inspect_parameter

    @inspect_parameter.setter
    def inspect_parameter(self, value: str | None) -> None:
        self._inspect_parameter = value

    @property
    def plot_style(self) -> str | None:
        return self._plot_style

    @plot_style.setter
    def plot_style(self, value: str | None) -> None:
        self._plot_style = value

    def __setattr__(self, key: str, value: Any) -> None:
        # Check if the class is frozen, allowing the internal '_frozen'
        # flag itself to be set
        if getattr(self, '_frozen', False) and key != '_frozen':
            raise RuntimeError(f"Cannot modify property '{key}': "
                               f"Configuration is frozen for workers!")

        # handle -s or --set
        if key == "set_parameter":
            if value and '=' in value:
                # Split into parameter key and value string
                raw_key, raw_val = value.split('=', 1)
                param_key = raw_key.strip()
                val_str = raw_val.strip()

                if param_key not in self.MAIN_PARAM_FILE_MAP:
                    raise ValueError(
                        f"Param '{param_key}' passed via -s/--set is not a "
                        f"valid Starbug2 parameter. Please check spelling")

                property_name, target_type = (
                    self.MAIN_PARAM_FILE_MAP[param_key])

                # Convert empty arguments or strings to standard configurations
                if val_str == "" or val_str.lower() == "none":
                    super().__setattr__(property_name, None)
                elif target_type is bool:
                    super().__setattr__(property_name, bool(int(val_str)))
                else:
                    super().__setattr__(property_name, target_type(val_str))
            return

        super().__setattr__(key, value)

    def __str__(self) -> str:
        """Returns a formatted, categorised string representation of the
           configuration."""

        # Helper to format array shapes or lengths instead of full objects
        def _fmt(val: object) -> str:
            if isinstance(val, list):
                return f"list(len={len(val)})"
            if hasattr(val, "shape"):  # NumPy arrays / Astropy tables
                return f"array(shape={getattr(val, 'shape')})"
            return repr(val)

        return (
            f"<{self.__class__.__name__} (\n"
            f"  [Meta]\n"
            f"    frozen={self._frozen}, verbose_logs={self._verbose_logs}, "
            f"    n_cores={self._n_cores}\n"
            f"  [Main Actions]\n"
            f"    aperture_phot={self._do_aperture_photometry}, "
            f"    bgd_estimate={self._do_bgd_estimate}, "
            f"star_detection={self._do_star_detection}, "
            f"matching={self._do_matching}, "
            f"phot_routine={self._do_photometry_routine},"
            f"bgd_subtraction={self._do_bgd_subtraction}\n"
            f"  [Files & Paths]\n"
            f"    param_file={self._param_file!r}, ap_file={self._ap_file!r}, "
            f"bgd_file={self._background_file!r},"
            f" output_file={self._output_file!r}, "
            f"fits_images={_fmt(self._fits_images)}\n"
            f"  [Instrument & Filters]\n"
            f"    hdu_name={self._hdu_name!r}, filter={self._filter!r},"
            f" detector={self._detector_name!r}\n"
            f"  [Photometry & Geometry]\n"
            f"    fwhm={self._full_width_half_max}, "
            f"sigma_sky={self._sigma_sky}, "
            f"sigma_source={self._sigma_source}, "
            f"ap_radius={self._aperture_phot_radius}, "
            f"sky_annulus=[{self._sky_annulus_inner_radius}, "
            f"{self._sky_annulus_outer_radius}]\n"
            f"  [Artificial Stars (AST)]\n"
            f"    do_ast={self._do_artificial_star_test}, "
            f"tests_count={self._artificial_star_tests_count}, "
            f"stars_per_test={self._stars_per_artificial_test}, "
            f"ast_loader={_fmt(self._ast_loader)}\n"
            f")>"
        )

    def __repr__(self) -> str:
        return self.__str__()
