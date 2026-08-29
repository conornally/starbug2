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
import math
import os
from multiprocessing import shared_memory
from multiprocessing.shared_memory import SharedMemory
from typing import Final

import numpy as np
import pytest
from astropy.table import Table
from photutils.psf import ImagePSF

from starbug2.core.star_bug_config import StarBugMainConfig
from starbug2.core.starbug_main import StarbugBase
from tests.generic import TEST_PATH_STR, TEST_PSF_FITS
from starbug2.command_line_interfaces.ast import ast_main
from starbug2.constants import ExitStates, TableColumn
from tests import generic

# main ast run
c: np.ndarray = np.array([0, 0, 0], dtype=np.int64)
share_memory: SharedMemory = (
    shared_memory.SharedMemory(create=True, size=c.nbytes))
loading_buffer: np.ndarray = np.ndarray(
    c.shape, dtype=c.dtype, buffer=share_memory.buf)
TEST_FILTER_STRING: Final[str] = "-s FILTER=F444W"


def run(s):
    return ast_main(
        s.split() + [generic.TEST_IMAGE_FITS], share_memory, loading_buffer
    )


def update_config_for_fake_stars_into_blank(config: StarBugMainConfig) -> None:
    config.unfreeze()
    config.custom_filter = generic.TEST_CUSTOM_FILTER
    config.fits_images = [generic.TEST_BLANK]
    config.psf_file_override = TEST_PSF_FITS
    config.do_star_detection = False
    config.do_aperture_photometry = False
    config.do_artificial_star_test = True
    config.ast_load_psf = True
    config.ast_seed = 42
    config.ast_no_psf_phot = False
    config.ast_no_background = False
    config.test_magnitude_bright_limit = 20
    config.test_magnitude_faint_limit = 22
    config.sigma_sky = 4
    config.sigma_source = 10
    config.stars_per_artificial_test = 1
    config.artificial_star_tests_count = 1
    config.ast_save_added_image = True
    config.ast_save_added_image_path = TEST_PATH_STR
    config.output_file = TEST_PATH_STR
    config.freeze()


def test_run_basic():
    generic.verify_test_data_exists()
    generic.clean()
    # run single core
    assert (run(
        f"starbug2-ast -N10 -S10 -sPSF_FILE={TEST_PSF_FITS} "
        f"--output={generic.TEST_PATH} {TEST_FILTER_STRING}") ==
            ExitStates.EXIT_SUCCESS)

    # run multi-core
    assert (run(
        f"starbug2-ast -N30 -S10 -n3 -sPSF_FILE={TEST_PSF_FITS} "
        f"--output={generic.TEST_PATH} {TEST_FILTER_STRING}") ==
            ExitStates.EXIT_SUCCESS)

    # run multi-core in /tmp
    assert run(
        f"starbug2-ast -N30 -S10 -n3 -sPSF_FILE={TEST_PSF_FITS} -o /tmp/"
        f" {TEST_FILTER_STRING}") == ExitStates.EXIT_SUCCESS
    generic.clean()


@pytest.mark.skipif(
    os.getenv("RUN_STAR_BUG_PRODUCTION_TESTS") is None or
    os.getenv("RUN_STAR_BUG_PRODUCTION_TESTS") == "false",
    reason="Harsh stress test locked out of normal development runs due to "
           "length of time to run, CPU resources required which nearly slags"
           " the machine."
)
def test_run_harsh_inputs():
    generic.clean()
    assert (run(
        f"starbug2-ast -N1 -S1000 {TEST_FILTER_STRING}") ==
            ExitStates.EXIT_SUCCESS)
    assert (run(
        f"starbug2-ast -N1000 -S1 {TEST_FILTER_STRING}") ==
            ExitStates.EXIT_SUCCESS)
    assert (run(
        f"starbug2-ast -N10 -S10 -n100 {TEST_FILTER_STRING}") ==
            ExitStates.EXIT_SUCCESS)
    assert run(
        f"starbug2-ast -N1000 -S1000 -n1000"
        f" {TEST_FILTER_STRING}") == ExitStates.EXIT_SUCCESS
    generic.clean()


def test_add_stars_logic():
    generic.verify_test_data_exists()
    generic.clean()

    # create blank fits file.
    config: StarBugMainConfig = StarBugMainConfig()
    generic.create_blank_fits()

    # create config
    config.unfreeze()
    config.ast_add_stars = True
    config.fits_images = [generic.TEST_BLANK]
    config.ast_load_psf = True
    config.psf_file_override = TEST_PSF_FITS
    config.ast_seed = 42
    config.stars_per_artificial_test = 1
    config.artificial_star_tests_count = 1
    config.ast_save_added_image = True
    config.ast_save_added_image_path = TEST_PATH_STR
    config.custom_filter = generic.TEST_CUSTOM_FILTER
    config.output_file = TEST_PATH_STR
    config.sigma_sky = 4
    config.sigma_source = 10
    config.freeze()

    entrance: StarbugBase = StarbugBase(
        config=config, f_name=generic.TEST_BLANK, ap_file=None, bkg_file=None)

    # extract image data before adding stars
    original_image_data: np.ndarray = entrance.main_image().data.copy()

    # execute add stars
    entrance.run_starbug()

    # extract locations and new image data.
    locations: Table | None = entrance.ast_star_source_list
    assert locations is not None
    image_data: np.ndarray = entrance.main_image().data.copy()

    # set to 0 all the values in both images where these fake star has affected
    # the image.
    x_coord: int = int(locations[0][TableColumn.X_0])
    y_coord: int = int(locations[0][TableColumn.Y_0])
    half_width_x: int = math.ceil(ImagePSF(entrance.psf).origin[0])
    half_width_y: int = math.ceil(ImagePSF(entrance.psf).origin[1])

    # Slice the mask at the star location and flip it to False (
    # do not check these pixels)
    comparison_mask = np.ones_like(original_image_data, dtype=bool)
    comparison_mask[
        y_coord - half_width_y: y_coord + half_width_y,
        x_coord - half_width_x: x_coord + half_width_x
    ] = False

    # test outside fake star location it is identical
    np.testing.assert_allclose(
        original_image_data[comparison_mask],
        image_data[comparison_mask], rtol=1e-5, atol=1e-4
    )

    # check inside the star it's not identical
    star_box_original = original_image_data[~comparison_mask]
    star_box_new = image_data[~comparison_mask]
    assert not np.array_equal(star_box_original, star_box_new)


def test_results_for_execute_as_test():
    generic.verify_test_data_exists()
    generic.clean()

    # create blank fits file.
    config: StarBugMainConfig = StarBugMainConfig()
    generic.create_blank_fits()
    config.unfreeze()
    config.do_star_detection = True
    config.do_aperture_photometry = True
    config.ast_no_psf_phot = False
    config.ast_no_background = False
    config.custom_filter = generic.TEST_CUSTOM_FILTER
    config.fits_images = [generic.TEST_BLANK]
    config.psf_file_override = TEST_PSF_FITS
    config.output_file = TEST_PATH_STR
    config.sigma_sky = 4
    config.sigma_source = 10
    config.freeze()

    entrance: StarbugBase = StarbugBase(
        config=config, f_name=generic.TEST_BLANK, ap_file=None, bkg_file=None)

    # execute add stars
    entrance.run_starbug()

    detections: Table | None = entrance.detections
    assert detections is not None
    background_detections: Table = detections.copy()

    # create config
    update_config_for_fake_stars_into_blank(config)

    entrance: StarbugBase = StarbugBase(
        config=config, f_name=generic.TEST_BLANK, ap_file=None, bkg_file=None)

    # execute add stars
    entrance.run_starbug()

    # get data for comparison
    artificial_stars_detections: Table | None = entrance.ast_detections
    fake_star_locations: Table | None = entrance.ast_star_source_list

    assert artificial_stars_detections is not None
    assert fake_star_locations is not None

    assert len(artificial_stars_detections) == 1
    assert len(fake_star_locations) == 1

    # compare detections.
    matching_pixel_tolerance = 1
    for row in artificial_stars_detections:
        orig_x = row[TableColumn.X_DET]
        orig_y = row[TableColumn.Y_DET]

        # Calculate the Euclidean distance from this original star
        # to all detections in the new artificial image
        distances = np.sqrt(
            (fake_star_locations[
                 TableColumn.X_0] - orig_x) ** 2 +
            (fake_star_locations[
                 TableColumn.Y_0] - orig_y) ** 2
        )

        # Verify that at least one detection lies within our tolerance limit
        match_exists = np.any(distances <= matching_pixel_tolerance)
        assert match_exists

        # check that the background stars were not part of the found star
        distances = np.sqrt(
            (background_detections[
                 TableColumn.X_CENTROID] - orig_x) ** 2 +
            (background_detections[
                 TableColumn.Y_CENTROID] - orig_y) ** 2
        )
        match_exists = np.any(distances <= matching_pixel_tolerance)
        assert not match_exists

    generic.clean()


def test_ast_output_data():
    generic.verify_test_data_exists()
    generic.clean()

    # create blank fits file.
    config: StarBugMainConfig = StarBugMainConfig()
    generic.create_blank_fits(max_value=1.0)
    update_config_for_fake_stars_into_blank(config)

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
    config.custom_filter = generic.TEST_CUSTOM_FILTER
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

    # clean setup
    generic.clean()


if __name__ == "__main__":
    # This allows you to run the harsh test directly.
    test_run_harsh_inputs()
