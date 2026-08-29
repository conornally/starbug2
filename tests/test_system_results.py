from multiprocessing import shared_memory
from multiprocessing.shared_memory import SharedMemory
from typing import Final

import numpy as np

from starbug2.command_line_interfaces.main import starbug_internal_main
from starbug2.constants import ExitStates
from starbug2.core.star_bug_config import StarBugMainConfig
from tests import generic
from tests.generic import (
    TEST_PATH, TEST_BLANK, TEST_PATH_STR, TEST_PSF_FITS, create_default_config)
import os


TEST_NGC_FITS: Final[str] = str(
    os.path.join(str(TEST_PATH), "ngc6822_F770W_i2d.fits"))


class TestSystemResults:

    @staticmethod
    def _assert_results(
            lines: list[str],
            expected_plain: int | None = None,
            expected_background: int | None = None,
            expected_con_vol: int | None = None,
            expected_cleaning: int | None = None,
            expected_total: int | None = None, ratio_low: float = 1,
            ratio_high: float = 1):
        """
        processes output and determines if it passes expectations.

        :param expected_plain: how many expected to find in plain search.
        :param expected_background: how many expected to find in background
                                    search.
        :param expected_con_vol: how many expected to find in con search.
        :param expected_cleaning: how many expected to be cleaned away.
        :param expected_total: how many expected to find in total.
        :param lines: the output lines.
        :param ratio_low: the ratio for acceptance in low mode.
        :param ratio_high: the ratio for acceptance in high mode.
        :return: None
        """
        plain_line: str | None = None
        background_line: str | None = None
        con_vol_line: str | None = None
        cleaning: str | None = None
        total_line: str | None = None

        for line in lines:
            if "[PLAIN] pass:" in line:
                plain_line = line
            if "[BGD2D] pass:" in line:
                background_line = line
            if "[CONVL] pass" in line:
                con_vol_line = line
            if "cleaning" in line:
                cleaning = line
            if "Total" in line:
                total_line = line

        assert (
            plain_line is not None and background_line is not None
            and con_vol_line is not None and cleaning is not None
            and total_line is not None)

        count_part: str = plain_line.split("pass: ")[1]
        plain_count: int = int(count_part.split()[0])

        count_part = background_line.split("pass: ")[1]
        background_count: int = int(count_part.split()[0])

        count_part = con_vol_line.split("pass: ")[1]
        con_vol_count: int = int(count_part.split()[0])

        count_part = cleaning.split("cleaning ")[1]
        cleaning_count: int = int(count_part.split()[0])

        count_part = total_line.split("Total: ")[1]
        total_line_count: int = int(count_part.split()[0])

        if expected_plain is not None:
            assert ((expected_plain * ratio_low) <= plain_count
                    <= (expected_plain * ratio_high))
        if expected_background is not None:
            assert ((expected_background * ratio_low) <= background_count
                    <= (expected_background * ratio_high))
        if expected_con_vol is not None:
            assert ((expected_con_vol * ratio_low) <= con_vol_count
                    <= (expected_con_vol * ratio_high))
        if expected_cleaning is not None:
            assert ((expected_cleaning * ratio_low) <= cleaning_count
                    <= (expected_cleaning * ratio_high))
        if expected_total is not None:
            assert ((expected_total * ratio_low) <= total_line_count
                    <= (expected_total * ratio_high))

    def test_detection_on_proper_fits_file(self, capsys):
        config: StarBugMainConfig = create_default_config()
        config.unfreeze()
        config.custom_filter = 'F770W'
        config.fits_images = [TEST_NGC_FITS]
        config.psf_file_override = TEST_PSF_FITS
        config.do_star_detection = True
        config.verbose_logs = True
        config.output_file = TEST_PATH_STR
        config.freeze()
        exit_state: int = starbug_internal_main(config)
        assert exit_state == ExitStates.EXIT_SUCCESS
        generic.clean()

        captured = capsys.readouterr()
        lines = captured.out.splitlines()
        self._assert_results(lines, 4461, 4508, 21854, 3433, 18421, 0.8, 1.2)

    def test_detection_on_artificial_stars(self, capsys):
        generic.clean()
        c: np.ndarray = np.array([0, 0, 0], dtype=np.int64)
        share_memory: SharedMemory = (
            shared_memory.SharedMemory(create=True, size=c.nbytes))
        loading_buffer: np.ndarray = np.ndarray(
            c.shape, dtype=c.dtype, buffer=share_memory.buf)

        # set up config for creating psf
        generic.make_psf_for_blank()

        # set up config for artificial stars
        config: StarBugMainConfig = create_default_config()
        config.unfreeze()
        config.custom_filter = 'F770W'
        config.fits_images = [TEST_BLANK]
        config.psf_file_override = TEST_PSF_FITS
        config.verbose_logs = True

        # config to set off detection with psf
        config.stars_per_artificial_test = 15
        config.ast_save_added_image = True
        config.ast_save_added_image_path = TEST_PATH_STR
        config.output_file = TEST_PATH_STR
        config.ast_add_stars = True
        config.ast_loader = loading_buffer
        config.ast_test_index = 0
        config.artificial_star_tests_count = 1
        config.ast_auto_save = 10
        config.ast_load_psf = True
        config.do_star_detection = True

        config.zero_point_magnitude = 25
        config.test_magnitude_bright_limit = 12
        config.test_magnitude_faint_limit = 13
        config.full_width_half_max = 5.0
        config.ricker_wavelet_radius = 2.5
        config.sigma_source = 100.0

        config.sigma_sky = 1.0
        # Ke
        config.aperture_phot_radius = 5.0
        config.sharp_cutoff_low = 0.4
        config.sharp_cutoff_high = 1.0
        config.round1_cutoff_high = 0.5
        config.round2_cutoff_high = 0.5

        # --- FIXED SOURCE EXTRACTION PARAMS ---
        config.do_convolution = True
        config.clean_sources = True
        config.sharp_cutoff_low = 0.2
        config.freeze()

        # create empty fits file
        generic.create_blank_fits(use_noise=False)

        # execute detection.
        exit_state: int = starbug_internal_main(config)
        assert exit_state == ExitStates.EXIT_SUCCESS

        # check results.
        captured = capsys.readouterr()
        lines = captured.out.splitlines()
        self._assert_results(
            lines, expected_total=15, ratio_high=2, ratio_low=-1)
        generic.clean()

    def test_artificial_star_residual(self, capsys) -> None:
        generic.clean()
        c: np.ndarray = np.array([0, 0, 0], dtype=np.int64)
        share_memory: SharedMemory = (
            shared_memory.SharedMemory(create=True, size=c.nbytes))
        loading_buffer: np.ndarray = np.ndarray(
            c.shape, dtype=c.dtype, buffer=share_memory.buf)

        # set up config for creating psf
        generic.make_psf_for_blank()

        # create empty fits file
        generic.create_blank_fits()

        # set up config for artificial stars
        config: StarBugMainConfig = create_default_config()
        config.unfreeze()
        config.custom_filter = generic.TEST_CUSTOM_FILTER
        config.fits_images = [TEST_BLANK]
        config.psf_file_override = TEST_PSF_FITS
        config.verbose_logs = True
        config.output_file = TEST_PATH_STR

        # config to set off detection with psf
        config.stars_per_artificial_test = 15
        config.ast_save_added_image = True
        config.ast_save_added_image_path = TEST_PATH_STR
        config.ast_add_stars = True
        config.ast_loader = loading_buffer
        config.do_star_detection = True
        config.ast_seed = 42
        config.test_magnitude_bright_limit = 20
        config.test_magnitude_faint_limit = 22
        config.sigma_sky = 4
        config.sigma_source = 10
        config.freeze()

        # execute detection / phot / residual.
        exit_state: int = starbug_internal_main(config)
        assert exit_state == ExitStates.EXIT_SUCCESS

        # check results.
        captured = capsys.readouterr()
        lines = captured.out.splitlines()
        self._assert_results(
            lines, expected_total=15, ratio_low=-1, ratio_high=2)
        generic.clean()
