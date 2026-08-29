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

from astropy.table import Table
from numpy.testing import assert_array_equal

from starbug2.constants import TableColumn
from starbug2.core.star_bug_config import StarBugMainConfig
from starbug2.core.starbug_main import StarbugBase
from tests.test_ast_magnitudes import FIXED_FILTER_FOR_SCIENCE_CONFIDENCE
from tests.generic import TEST_PATH_STR, TEST_PSF_FITS
from tests import generic
import numpy as np


def update_config_for_fake_stars_into_image_fits(
        config: StarBugMainConfig, use_psf: bool) -> None:
    """
    this populates a config with the default values needed to run artificial
    stars as we want
    :param config: the config to adjust
    :type config: StarBugMainConfig
    :param use_psf: use psf or not
    :type use_psf: bool
    :return: None
    :rtype: None
    """
    config.unfreeze()
    config.custom_filter = FIXED_FILTER_FOR_SCIENCE_CONFIDENCE
    config.fits_images = [generic.TEST_IMAGE_FITS]
    config.psf_file_override = TEST_PSF_FITS
    config.do_star_detection = False
    config.do_aperture_photometry = False
    config.do_artificial_star_test = True
    config.ast_load_psf = True
    config.ast_seed = 42
    config.sigma_sky = 4
    config.sigma_source = 10
    config.test_magnitude_bright_limit = 20
    config.test_magnitude_faint_limit = 22
    config.stars_per_artificial_test = 1
    config.artificial_star_tests_count = 1
    config.ast_no_background = use_psf
    config.ast_no_psf_phot = use_psf
    config.ast_save_added_image = True
    config.ast_save_added_image_path = TEST_PATH_STR
    config.output_file = TEST_PATH_STR
    config.do_artificial_star_test_results = True
    config.plot_ast = os.path.join(TEST_PATH_STR, "plot")
    config.ast_plot_filename = os.path.join(TEST_PATH_STR, "plot")
    config.ast_save_added_image_path = TEST_PATH_STR
    config.freeze()


def update_config_for_basic_into_image_fits(
        config: StarBugMainConfig, use_psf: bool) -> None:
    """
    this populates a config with the default values needed to run basic
    detect and aperture
    :param config: the config to adjust
    :type config: StarBugMainConfig
    :param use_psf: use psf or not
    :type use_psf: bool
    :return: None
    :rtype: None
    """
    config.unfreeze()
    config.custom_filter = FIXED_FILTER_FOR_SCIENCE_CONFIDENCE
    config.fits_images = [
        os.path.join(TEST_PATH_STR, "inserted_image_for_test_0.fits")]
    config.psf_file_override = TEST_PSF_FITS
    config.do_star_detection = True
    config.do_aperture_photometry = True
    config.ast_load_psf = True
    config.sigma_sky = 4
    config.sigma_source = 10
    config.do_bgd_estimate = True
    config.do_photometry_routine = use_psf
    config.test_magnitude_bright_limit = 20
    config.test_magnitude_faint_limit = 22
    config.output_file = TEST_PATH_STR
    config.freeze()


def find_mask_for_artificial_star(
        fake_star_locations: Table, ap_file: Table, x_label: str) -> int:
    """
    figures which row to filter out of the main detections to acquire the
    fake star data.
    :param fake_star_locations: the fake star locations
    :type fake_star_locations: Table.
    :param ap_file: the detections table.
    :type ap_file: Table
    :param x_label: the label for the ap file column
    :type x_label: str
    :return: index that identifies the correct fake star location.
    :rtype: int
    """
    target_x = fake_star_locations[TableColumn.X_0][0]
    return int(np.argmin(np.abs(ap_file[x_label] - target_x)))


def test_artificial_stars_vs_basic():
    """
    tests that a run with artificial stars behaves the same as a basic run.

    For example. once stars are added. the results from detection and aperture
    produce the same as the artificial star run with the same image.
    :return: None
    """
    generic.verify_test_data_exists()
    generic.clean()

    config: StarBugMainConfig = StarBugMainConfig()
    update_config_for_fake_stars_into_image_fits(config, True)

    entrance: StarbugBase = StarbugBase(
        config=config, f_name=generic.TEST_IMAGE_FITS, ap_file=None,
        bkg_file=None)
    # execute add stars and do artificial test
    entrance.run_starbug()

    artificial_stars_detections: Table | None = entrance.ast_detections
    fake_star_locations: Table | None = entrance.ast_star_source_list
    output_file: str = os.path.join(TEST_PATH_STR, "image-ast.fits")

    assert os.path.exists(output_file)
    ast_file: Table = Table.read(
        os.path.join(TEST_PATH_STR, "image-ast.fits"),
        format="fits", hdu=2).copy()
    ap_file: Table = Table.read(
        os.path.join(TEST_PATH_STR, "image-ap.fits"), hdu=1).copy()

    # execute basic detection and aperture.
    config: StarBugMainConfig = StarBugMainConfig()
    update_config_for_basic_into_image_fits(config, False)
    entrance: StarbugBase = StarbugBase(
        config=config,
        f_name=os.path.join(TEST_PATH_STR, "inserted_image_for_test_0.fits"),
        ap_file=None,
        bkg_file=None)
    entrance.run_starbug()

    basic_stars_detections: Table | None = entrance.detections
    basic_ap_file: Table = Table.read(
        os.path.join(TEST_PATH_STR, "inserted_image_for_test_0-ap.fits"),
        hdu=1).copy()

    # assert the psf files do not exist
    assert (not os.path.exists(
        os.path.join(TEST_PATH_STR, "inserted_image_for_test_0.psf")))

    # check table contents
    assert (ap_file.colnames == basic_ap_file.colnames)
    assert_array_equal(ap_file[TableColumn.X_CENTROID],
                       basic_ap_file[TableColumn.X_CENTROID])
    assert_array_equal(ap_file[TableColumn.RA], basic_ap_file[TableColumn.RA])
    assert_array_equal(
        ap_file[TableColumn.DEC], basic_ap_file[TableColumn.DEC])
    assert_array_equal(
        ap_file[TableColumn.SHARPNESS], basic_ap_file[TableColumn.SHARPNESS])
    assert_array_equal(
        ap_file[TableColumn.ROUNDNESS1], basic_ap_file[TableColumn.ROUNDNESS1])
    assert_array_equal(
        ap_file[TableColumn.ROUNDNESS2], basic_ap_file[TableColumn.ROUNDNESS2])
    assert_array_equal(
        ap_file[TableColumn.SMOOTHNESS], basic_ap_file[TableColumn.SMOOTHNESS])
    assert_array_equal(
        ap_file[TableColumn.SKY], basic_ap_file[TableColumn.SKY])
    assert_array_equal(
        ap_file[TableColumn.Y_CENTROID], basic_ap_file[TableColumn.Y_CENTROID])
    assert_array_equal(
        ap_file[TableColumn.FLUX],  basic_ap_file[TableColumn.FLUX])
    assert_array_equal(
        ap_file[TableColumn.E_FLUX], basic_ap_file[TableColumn.E_FLUX])
    assert_array_equal(
        ap_file[FIXED_FILTER_FOR_SCIENCE_CONFIDENCE],
        basic_ap_file[FIXED_FILTER_FOR_SCIENCE_CONFIDENCE])
    assert_array_equal(
        ap_file[f"e{FIXED_FILTER_FOR_SCIENCE_CONFIDENCE}"],
        basic_ap_file[f"e{FIXED_FILTER_FOR_SCIENCE_CONFIDENCE}"])

    # ensure data exists
    assert artificial_stars_detections is not None
    assert fake_star_locations is not None
    assert basic_stars_detections is not None

    # verify
    assert len(artificial_stars_detections) == 1
    assert len(fake_star_locations) == 1
    assert len(basic_stars_detections) == 10

    # find location
    fake_star_index: int = find_mask_for_artificial_star(
        fake_star_locations, ap_file, TableColumn.X_CENTROID)

    # verify mag from ast
    assert_array_equal(
        ast_file[TableColumn.MAG_DET],
        basic_ap_file[FIXED_FILTER_FOR_SCIENCE_CONFIDENCE][fake_star_index]
    )
    assert_array_equal(
        ast_file[TableColumn.MAG_DET],
        ap_file[FIXED_FILTER_FOR_SCIENCE_CONFIDENCE][fake_star_index]
    )
    assert_array_equal(
        ast_file[TableColumn.X_DET],
        basic_ap_file[TableColumn.X_CENTROID][fake_star_index]
    )
    assert_array_equal(
        ast_file[TableColumn.Y_DET],
        basic_ap_file[TableColumn.Y_CENTROID][fake_star_index])

    generic.clean()


def test_artificial_stars_vs_basic_psf():
    """
    tests that a run with artificial stars behaves the same as a basic run with
    psf photometry.

    For example. once stars are added. the results from detection, aperture,
    psf photometry produce the same as the artificial star run with the same
    image.
    :return: None
    """
    generic.verify_test_data_exists()
    generic.clean()

    config: StarBugMainConfig = StarBugMainConfig()
    update_config_for_fake_stars_into_image_fits(config, False)

    entrance: StarbugBase = StarbugBase(
        config=config, f_name=generic.TEST_IMAGE_FITS, ap_file=None,
        bkg_file=None)
    # execute add stars and do artificial test
    entrance.run_starbug()

    artificial_stars_detections: Table | None = entrance.ast_detections
    ast_basic_detections: Table | None = entrance.detections
    fake_star_locations: Table | None = entrance.ast_star_source_list
    output_file: str = os.path.join(TEST_PATH_STR, "image-ast.fits")

    # find location
    assert ast_basic_detections is not None
    assert fake_star_locations is not None
    fake_star_index = find_mask_for_artificial_star(
        fake_star_locations, ast_basic_detections,
        TableColumn.X_CENTROID)

    assert os.path.exists(output_file)
    ast_ast_file: Table = Table.read(
        os.path.join(TEST_PATH_STR, "image-ast.fits"),
        format="fits", hdu=2).copy()
    ast_ap_file: Table = Table.read(
        os.path.join(TEST_PATH_STR, "image-ap.fits"), hdu=1).copy()
    ast_psf_file: Table = Table.read(
        os.path.join(TEST_PATH_STR, "image-psf.fits"), hdu=1).copy()

    # execute basic detection and aperture.
    config: StarBugMainConfig = StarBugMainConfig()
    update_config_for_basic_into_image_fits(config, True)
    entrance: StarbugBase = StarbugBase(
        config=config,
        f_name=os.path.join(TEST_PATH_STR, "inserted_image_for_test_0.fits"),
        ap_file=None,
        bkg_file=None)
    entrance.run_starbug()

    basic_stars_detections: Table | None = entrance.detections
    basic_ap_file: Table = Table.read(
        os.path.join(TEST_PATH_STR, "inserted_image_for_test_0-ap.fits"),
        hdu=1).copy()
    basic_psf_file: Table = Table.read(
        os.path.join(TEST_PATH_STR, "inserted_image_for_test_0-psf.fits"),
        hdu=1).copy()

    # check ap tables contents
    assert (ast_ap_file.colnames == basic_ap_file.colnames)
    assert_array_equal(ast_ap_file[TableColumn.X_CENTROID],
                       basic_ap_file[TableColumn.X_CENTROID])
    assert_array_equal(
        ast_ap_file[TableColumn.RA], basic_ap_file[TableColumn.RA])
    assert_array_equal(
        ast_ap_file[TableColumn.DEC], basic_ap_file[TableColumn.DEC])
    assert_array_equal(
        ast_ap_file[TableColumn.SHARPNESS],
        basic_ap_file[TableColumn.SHARPNESS])
    assert_array_equal(
        ast_ap_file[TableColumn.ROUNDNESS1],
        basic_ap_file[TableColumn.ROUNDNESS1])
    assert_array_equal(
        ast_ap_file[TableColumn.ROUNDNESS2],
        basic_ap_file[TableColumn.ROUNDNESS2])
    assert_array_equal(
        ast_ap_file[TableColumn.SMOOTHNESS],
        basic_ap_file[TableColumn.SMOOTHNESS])
    assert_array_equal(
        ast_ap_file[TableColumn.SKY], basic_ap_file[TableColumn.SKY])
    assert_array_equal(
        ast_ap_file[TableColumn.Y_CENTROID],
        basic_ap_file[TableColumn.Y_CENTROID])
    assert_array_equal(
        ast_ap_file[TableColumn.FLUX],  basic_ap_file[TableColumn.FLUX])
    assert_array_equal(
        ast_ap_file[TableColumn.E_FLUX], basic_ap_file[TableColumn.E_FLUX])
    assert_array_equal(
        ast_ap_file[FIXED_FILTER_FOR_SCIENCE_CONFIDENCE],
        basic_ap_file[FIXED_FILTER_FOR_SCIENCE_CONFIDENCE])
    assert_array_equal(
        ast_ap_file[f"e{FIXED_FILTER_FOR_SCIENCE_CONFIDENCE}"],
        basic_ap_file[f"e{FIXED_FILTER_FOR_SCIENCE_CONFIDENCE}"])

    # ensure data exists
    assert artificial_stars_detections is not None
    assert fake_star_locations is not None
    assert basic_stars_detections is not None
    assert ast_basic_detections is not None

    # verify table sizes
    assert len(artificial_stars_detections) == 1
    assert len(fake_star_locations) == 1
    assert len(basic_stars_detections) == 10
    assert len(ast_basic_detections) == 10

    # verify mag from ast
    assert_array_equal(
        ast_ast_file[TableColumn.MAG_DET],
        basic_psf_file[FIXED_FILTER_FOR_SCIENCE_CONFIDENCE][fake_star_index]
    )
    assert_array_equal(
        ast_ast_file[TableColumn.MAG_DET],
        ast_psf_file[FIXED_FILTER_FOR_SCIENCE_CONFIDENCE][fake_star_index]
    )
    assert_array_equal(
        ast_ast_file[TableColumn.X_DET],
        basic_psf_file[TableColumn.X_INIT][fake_star_index]
    )
    assert_array_equal(
        ast_ast_file[TableColumn.Y_DET],
        basic_psf_file[TableColumn.Y_INIT][fake_star_index])

    # verify psf files are identical.
    assert_array_equal(
        ast_psf_file[TableColumn.X_INIT], basic_psf_file[TableColumn.X_INIT])
    assert_array_equal(
        ast_psf_file[TableColumn.Y_INIT], basic_psf_file[TableColumn.Y_INIT])
    assert_array_equal(
        ast_psf_file[TableColumn.RA], basic_psf_file[TableColumn.RA])
    assert_array_equal(
        ast_psf_file[TableColumn.DEC], basic_psf_file[TableColumn.DEC])
    assert_array_equal(
        ast_psf_file[TableColumn.FLAG], basic_psf_file[TableColumn.FLAG])
    assert_array_equal(
        ast_psf_file[f"ap_{FIXED_FILTER_FOR_SCIENCE_CONFIDENCE}"],
        basic_psf_file[f"ap_{FIXED_FILTER_FOR_SCIENCE_CONFIDENCE}"])
    assert_array_equal(
        ast_psf_file[TableColumn.X_FIT], basic_psf_file[TableColumn.X_FIT])
    assert_array_equal(
        ast_psf_file[TableColumn.Y_FIT], basic_psf_file[TableColumn.Y_FIT])
    assert_array_equal(
        ast_psf_file[TableColumn.FLUX], basic_psf_file[TableColumn.FLUX])
    assert_array_equal(
        ast_psf_file[TableColumn.E_FLUX], basic_psf_file[TableColumn.E_FLUX])
    assert_array_equal(
        ast_psf_file[TableColumn.XY_DEV], basic_psf_file[TableColumn.XY_DEV])
    assert_array_equal(
        ast_psf_file[TableColumn.Q_FIT], basic_psf_file[TableColumn.Q_FIT])
    assert_array_equal(
        ast_psf_file[FIXED_FILTER_FOR_SCIENCE_CONFIDENCE],
        basic_psf_file[FIXED_FILTER_FOR_SCIENCE_CONFIDENCE])
    assert_array_equal(
        ast_psf_file[f"e{FIXED_FILTER_FOR_SCIENCE_CONFIDENCE}"],
        basic_psf_file[f"e{FIXED_FILTER_FOR_SCIENCE_CONFIDENCE}"])

    generic.clean()
