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
from typing import Final

import numpy as np
import pytest
from astropy.table import Table

from starbug2.constants import ExitStates, TableColumn
from starbug2.core.main_components.one_time_runs import starbug_one_time_runs
from starbug2.core.star_bug_config import StarBugMainConfig
from starbug2.core.starbug_main import StarbugBase
from tests.generic import TEST_PATH_STR, TEST_PSF_FITS
from tests import generic

FIXED_FILTER_FOR_SCIENCE_CONFIDENCE = "F444W"

# these values come from an execution of starbug2 on a MacBook, as grounding
# values for verifying that the newest starbug 2 generates the same values.
X_COORD_LOCATION_OLD_STARBUG: Final = 77.74946
DETECTED_MAG_OLD_STARBUG: Final = 21.78636
DETECTED_ERROR_OLD_STARBUG: Final = 0.09465


def update_config_for_fake_stars_into_image_fits(
        config: StarBugMainConfig,
        custom_filter: str = FIXED_FILTER_FOR_SCIENCE_CONFIDENCE) -> None:
    """
    this populates a config with the default values needed to run artificial
    stars as we want
    :param config: the config to adjust
    :type config: StarBugMainConfig
    :param custom_filter: the custom filter
    :type custom_filter: str
    :return: None
    :rtype: None
    """
    config.unfreeze()
    config.custom_filter = custom_filter
    config.fits_images = [generic.TEST_IMAGE_FITS]
    config.psf_file_override = TEST_PSF_FITS
    config.do_star_detection = False
    config.do_aperture_photometry = False
    config.do_artificial_star_test = True
    config.ast_load_psf = True
    config.ast_seed = 42
    config.test_magnitude_bright_limit = 20
    config.test_magnitude_faint_limit = 22
    config.stars_per_artificial_test = 1
    config.artificial_star_tests_count = 1
    config.ast_save_added_image = True
    config.ast_save_added_image_path = TEST_PATH_STR
    config.output_file = TEST_PATH_STR
    config.freeze()


def test_ast_output_data():
    """
    tests that artificial stars can generate the correct output utilising
    a test fits file.
    :return: None
    """
    generic.verify_test_data_exists()
    generic.clean()

    # create blank fits file.
    config: StarBugMainConfig = StarBugMainConfig()
    generic.create_blank_fits()
    update_config_for_fake_stars_into_image_fits(
        config, generic.TEST_CUSTOM_FILTER)
    config.unfreeze()
    config.verbose_logs = True
    config.sigma_sky = 4
    config.sigma_source = 10
    config.generate_residual_image = True
    config.ast_no_psf_phot = False
    config.ast_no_background = False
    config.freeze()

    entrance: StarbugBase = StarbugBase(
        config=config, f_name=generic.TEST_BLANK, ap_file=None, bkg_file=None)

    # execute add stars and do test
    entrance.run_starbug()

    artificial_stars_detections: Table | None = entrance.ast_detections
    fake_star_locations: Table | None = entrance.ast_star_source_list

    assert artificial_stars_detections is not None
    assert fake_star_locations is not None
    assert len(artificial_stars_detections) == 1
    assert len(fake_star_locations) == 1

    assert (artificial_stars_detections[0][TableColumn.X_DET] ==
            pytest.approx(fake_star_locations[0][TableColumn.X_0], abs=0.5))
    assert (artificial_stars_detections[0][TableColumn.Y_DET] ==
            pytest.approx(fake_star_locations[0][TableColumn.Y_0], abs=0.5))
    assert (artificial_stars_detections[0][TableColumn.FLUX_DET] ==
            pytest.approx(fake_star_locations[0][TableColumn.FLUX], abs=0.1))

    # execute output generation
    config = StarBugMainConfig()
    config.custom_filter = FIXED_FILTER_FOR_SCIENCE_CONFIDENCE
    config.fits_images = [generic.TEST_BLANK]
    config.psf_file_override = TEST_PSF_FITS
    config.do_artificial_star_test_results = True
    config.plot_ast = os.path.join(TEST_PATH_STR, "plot")
    config.ast_plot_filename = os.path.join(TEST_PATH_STR, "plot")
    config.ast_out_tables = entrance.ast_test_results
    config.ast_save_added_image_path = TEST_PATH_STR
    config.output_file = TEST_PATH_STR

    entrance: StarbugBase = StarbugBase(
        config=config, f_name=generic.TEST_BLANK, ap_file=None, bkg_file=None)
    entrance.run_starbug()

    # check output generated.
    output_file: str = os.path.join(TEST_PATH_STR, "blank-ast.fits")
    output_file2: str = os.path.join(TEST_PATH_STR, "plot.png")
    assert os.path.exists(output_file)
    assert os.path.exists(output_file2)

    # generate local param file for sending back to starbug1
    config.generate_local_param_file = True
    starbug_one_time_runs(config)

    # clean setup
    generic.clean()


def test_ast_output_psf_photo_data():
    """
    this test uses the image.fits inside test folder and sees if it can
    generate values which match with old starbug.
    :return: None
    """
    generic.verify_test_data_exists()
    generic.clean()

    # create blank fits file.
    config: StarBugMainConfig = StarBugMainConfig()
    generic.create_blank_fits()
    update_config_for_fake_stars_into_image_fits(config)
    config.unfreeze()
    config.verbose_logs = True
    config.sigma_sky = 4
    config.sigma_source = 10
    config.ast_no_background = False
    config.ast_no_psf_phot = False
    config.generate_residual_image = True
    config.freeze()

    entrance: StarbugBase = StarbugBase(
        config=config, f_name=generic.TEST_IMAGE_FITS, ap_file=None,
        bkg_file=None)

    # execute add stars and do test
    entrance.run_starbug()

    artificial_stars_detections: Table | None = entrance.ast_detections
    fake_star_locations: Table | None = entrance.ast_star_source_list

    assert artificial_stars_detections is not None
    assert fake_star_locations is not None
    assert len(artificial_stars_detections) == 1
    assert len(fake_star_locations) == 1

    assert (artificial_stars_detections[0][TableColumn.X_DET] ==
            pytest.approx(fake_star_locations[0][TableColumn.X_0], abs=0.5))
    assert (artificial_stars_detections[0][TableColumn.Y_DET] ==
            pytest.approx(fake_star_locations[0][TableColumn.Y_0], abs=0.5))
    assert (artificial_stars_detections[0][TableColumn.FLUX_DET] ==
            pytest.approx(fake_star_locations[0][TableColumn.FLUX], abs=0.1))

    # execute output generation
    config = StarBugMainConfig()
    config.custom_filter = generic.TEST_CUSTOM_FILTER
    config.fits_images = [generic.TEST_IMAGE_FITS]
    config.psf_file_override = TEST_PSF_FITS
    config.do_artificial_star_test_results = True
    config.plot_ast = os.path.join(TEST_PATH_STR, "plot")
    config.ast_plot_filename = os.path.join(TEST_PATH_STR, "plot")
    config.ast_out_tables = entrance.ast_test_results
    config.ast_save_added_image_path = TEST_PATH_STR
    config.output_file = TEST_PATH_STR

    entrance: StarbugBase = StarbugBase(
        config=config, f_name=generic.TEST_IMAGE_FITS, ap_file=None,
        bkg_file=None)
    entrance.run_starbug()

    # check output generated.
    output_file: str = os.path.join(TEST_PATH_STR, "image-ast.fits")
    output_file2: str = os.path.join(TEST_PATH_STR, "plot.png")
    assert os.path.exists(output_file)
    assert os.path.exists(output_file2)

    # check param file can be generated and exists.
    config.generate_local_param_file = True
    exit_state: ExitStates = starbug_one_time_runs(config)
    assert (exit_state == ExitStates.EXIT_SUCCESS)
    param_file_name: str = os.path.join(TEST_PATH_STR, "starbug.param")
    assert (os.path.isfile(param_file_name))

    # check output values are sensible.
    ast_file: Table = Table.read(
        os.path.join(TEST_PATH_STR, "image-ast.fits"), format="fits", hdu=2)
    assert (ast_file[TableColumn.MAG_DET] ==
            pytest.approx(DETECTED_MAG_OLD_STARBUG, 0.1))

    ap_file: Table = Table.read(
        os.path.join(TEST_PATH_STR, "image-ap.fits"), hdu=1)
    x_centroids = ap_file[TableColumn.X_CENTROID]
    mask = np.isclose(x_centroids, X_COORD_LOCATION_OLD_STARBUG, atol=1e-5)
    matching_rows = ap_file[mask]
    assert len(matching_rows) == 1

    # NOTE numbers were determined from old starbug2 from connor
    assert (matching_rows[FIXED_FILTER_FOR_SCIENCE_CONFIDENCE] ==
            pytest.approx(DETECTED_MAG_OLD_STARBUG, abs=1e-5))
    assert (matching_rows[f"e{FIXED_FILTER_FOR_SCIENCE_CONFIDENCE}"] ==
            pytest.approx(DETECTED_ERROR_OLD_STARBUG, abs=1e-5))
    # clean setup
    generic.clean()
