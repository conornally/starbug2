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
import glob
from urllib import request
from typing import Final

import numpy as np
import pytest
from astropy.io import fits

from starbug2.command_line_interfaces.main import starbug_internal_main
from starbug2.constants import (
    STAR_BUG_TEST_DAT_ENV, ImageHeaderTags, MIRI_STRING, MIRI_IMAGE,
    DEFAULT_BUNIT, HeaderTags)
from starbug2.jwst_support.initialise_psf_data import download_ap_corr_files
from starbug2.core.star_bug_config import StarBugMainConfig
from starbug2.utilities.utils import get_data_path

# paths to test files
TEST_PATH: Final[str | None] = os.getenv(STAR_BUG_TEST_DAT_ENV)
if TEST_PATH is None:
    raise Exception("cant find the test data environmental variable")
TEST_PATH_STR: Final[str] = str(TEST_PATH)
TEST_IMAGE_FITS: Final[str] = os.path.join(TEST_PATH, "image.fits")
TEST_PSF_FITS: Final[str] = os.path.join(TEST_PATH, "psf.fits")
TEST_NGC_FITS: Final[str] = os.path.join(TEST_PATH, "ngc6822_F770W_i2d.fits")
TEST_README: Final[str] = os.path.join(TEST_PATH, "readme.txt")
TEST_BLANK: Final[str] = str(os.path.join(str(TEST_PATH), "blank.fits"))
TEST_AST_FILLED: Final[str] = str(
    os.path.join(str(TEST_PATH), "inserted_image_for_test_0.fits"))
TEST_SEED = 42

# the filter string for tests to ensure they all use the same stuff
TEST_CUSTOM_FILTER = "F770W"
TEST_FILTER_STRING_NO_G = "-s FILTER=F444W"
TEST_FILTER_STRING = "-s FILTER=F444W -G"
GITHUB_RELEASE_URL = (
    "https://github.com/alan-stokes/starbug2/releases/download/TEST_DATA/")
REQUIRED_FILES = [
    "image.fits", "psf.fits", f"ngc6822_{TEST_CUSTOM_FILTER}_i2d.fits"]


def verify_test_data_exists() -> None:
    # Check if the specific FITS file is missing
    if not (os.path.exists(TEST_IMAGE_FITS)
            and os.path.exists(TEST_PSF_FITS)
            and os.path.exists(TEST_NGC_FITS)):
        print(
            "\n⚠️ Test file missing due to merge. "
            "Downloading all from GitHub Releases...")

        for filename in REQUIRED_FILES:
            file_path = os.path.join(str(TEST_PATH), filename)
            url = f"{GITHUB_RELEASE_URL}/{filename}"
            if not os.path.exists(file_path):
                try:
                    request.urlretrieve(url, file_path)
                except Exception as e:
                    pytest.fail(
                        f"Failed to download test asset from GitHub "
                        f"Release: {e}")


def create_default_config() -> StarBugMainConfig:
    """
    creates a default config where everything points to the test output dir
    :return: a config
    :rtype StarBugMainConfig
    """
    config: StarBugMainConfig = StarBugMainConfig()
    config.output_file = TEST_PATH
    return config


def clean() -> None:
    """
    cleans up the test data folder for new tests.
    :return: None
    """
    files = glob.glob(os.path.join(str(TEST_PATH), "*"))
    files.remove(TEST_IMAGE_FITS)
    files.remove(TEST_PSF_FITS)
    files.remove(TEST_NGC_FITS)
    files.remove(TEST_README)
    for file_name in files:
        os.remove(file_name)
    if os.path.exists("dat/starbug.param"):
        os.remove("dat/starbug.param")


def check_shape(c, out) -> None:
    """
    checks shape.
    :param c: array 1
    :param out: array 2
    :return: None
    """
    assert np.shape(c) == np.shape(out)
    for m in range(len(c)):
        for n in range(len(c[m])):
            a = c[m][n]
            b = out[m][n]
            assert np.isnan(a) == np.isnan(b)
            if not np.isnan(a) or not np.isnan(b):
                assert a == b


def create_blank_fits(
        size: tuple[int, int] = (2048, 2048), use_noise: bool = True,
        max_value: float = 0.001):
    """
    creates a blank fits file.
    :param size: the size of the fits file.
    :type size: Tuple[int, int]
    :param use_noise: bool to test if we use background noise
    :type use_noise: false
    :param max_value: highest value in the background.
    :type max_value: float
    :return: None
    """
    print(f"Generating blank space image of size {size[0]}x{size[1]}...")

    # Create a 2D numpy array of zeros (using float32 for standard precision)
    blank_data = np.zeros(size, dtype=np.float64)

    # Create background noise: mean of 10.0 counts, standard deviation of 1.0
    rng = np.random.default_rng(seed=TEST_SEED)
    if use_noise:
        raw_noise = rng.normal(
            loc=max_value, scale=0.10, size=blank_data.shape)
        background_noise = np.clip(raw_noise, a_min=0.0, a_max=None)
    else:
        background_noise = np.zeros(size, dtype=np.float64)

    # Wrap the data inside a Primary HDU
    primary_hdu = fits.PrimaryHDU()

    # Add essential metadata headers so pipeline loaders don't choke
    header = primary_hdu.header
    header["EXTNAME"] = "PRIMARY"
    header["OBJECT"] = "BLANK_SPACE_CI"
    header["COMMENT"] = (
        "Artificial black space for starbug2 integration tests.")
    header[ImageHeaderTags.DETECTOR] = MIRI_IMAGE
    header[ImageHeaderTags.INSTRUMENT] = MIRI_STRING
    header[ImageHeaderTags.BUN_IT] = DEFAULT_BUNIT
    header[ImageHeaderTags.PIXAR_SR] = 9.31e-14
    header[ImageHeaderTags.PIXAR_A2] = 0.00396
    header[ImageHeaderTags.TELESCOPE] = "JWST"
    header[ImageHeaderTags.JWST] = 1
    header[HeaderTags.FILTER] = "F770W"

    # 2. SCI Extension (Science Data)
    sci_hdu = fits.ImageHDU(data=background_noise, name="SCI")
    # Copy relevant headers to SCI HDU if required by starbug
    sci_hdu.header[ImageHeaderTags.BUN_IT] = DEFAULT_BUNIT
    sci_hdu.header[ImageHeaderTags.DETECTOR] = MIRI_IMAGE
    sci_hdu.header[ImageHeaderTags.INSTRUMENT] = MIRI_STRING
    sci_hdu.header[ImageHeaderTags.PIXAR_SR] = 9.31e-14
    sci_hdu.header[ImageHeaderTags.PIXAR_A2] = 0.00396

    # 3. ERR Extension (Uncertainties / Error map)
    err_data = np.sqrt(np.abs(background_noise))
    err_hdu = fits.ImageHDU(data=err_data, name="ERR")

    # 4. DQ Extension (Data Quality flags - 0 means good pixel)
    dq_data = np.zeros(size, dtype=np.int32)
    dq_hdu = fits.ImageHDU(data=dq_data, name="DQ")  # type: ignore

    # 5. AREA Extension (Pixel area map - set to 1.0 or nominal pixel scale)
    area_data = np.ones(size, dtype=np.float32)
    area_hdu = fits.ImageHDU(data=area_data, name="AREA")

    # Combine all HDUs into an HDUList
    hdu_list: fits.HDUList = fits.HDUList(
        [primary_hdu, sci_hdu, err_hdu, dq_hdu, area_hdu])

    # Write the file out to disk
    # overwrite=True ensures test scripts can recreate this file on
    # every run
    hdu_list.writeto(TEST_BLANK, overwrite=True)
    print(f"Successfully saved to {TEST_BLANK}")


def make_psf_for_blank() -> None:
    """
    creates the psf used by blank.psf and downloads the ap_corr files
    :return: None
    """
    file_path: str = os.path.join(
        get_data_path(), f"{TEST_CUSTOM_FILTER}.fits")
    if os.path.exists(file_path):
        return

    psf_config: StarBugMainConfig = create_default_config()
    psf_config.custom_filter = TEST_CUSTOM_FILTER
    psf_config.generate_psf = True
    psf_config.detector_name = None
    psf_config.psf_fit_size = None
    starbug_internal_main(psf_config)
    download_ap_corr_files(get_data_path())
