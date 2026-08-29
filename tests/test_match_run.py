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

import pytest

from starbug2.command_line_interfaces.main import starbug_main
from starbug2.command_line_interfaces.match import match_main
from starbug2.constants import ExitStates
from tests.generic import (
    clean, TEST_IMAGE_FITS, TEST_FILTER_STRING, TEST_PATH_STR,
    verify_test_data_exists, TEST_PATH, TEST_PSF_FITS)

OUT_1_FITS: Final[str] = str(os.path.join(TEST_PATH_STR, "out1.fits"))
OUT_2_FITS: Final[str] = os.path.join(TEST_PATH_STR, "out2.fits")
OUT_1_AP_FITS: Final[str] = os.path.join(TEST_PATH_STR, "out1-ap.fits")
OUT_2_AP_FITS: Final[str] = os.path.join(TEST_PATH_STR, "out2-ap.fits")
IMAGE_AP_FITS: Final[str] = os.path.join(TEST_PATH_STR, "image-ap.fits")


def run(s):
    return match_main(s.split())


def test_match_start():
    verify_test_data_exists()
    assert (run(f"starbug2-match --output={TEST_PATH}") ==
            ExitStates.EXIT_EARLY)
    assert (run(f"starbug2-match -h --output={TEST_PATH}") ==
            ExitStates.EXIT_SUCCESS)
    assert (run(f"starbug2-match -vh --output={TEST_PATH}") ==
            ExitStates.EXIT_SUCCESS)


def test_match_bad_input():
    assert (run(f"starbug2-match --output={TEST_PATH}") ==
            ExitStates.EXIT_EARLY)
    assert (run(f"starbug2-match {TEST_IMAGE_FITS} --output={TEST_PATH}") ==
            ExitStates.EXIT_EARLY)
    assert (run(f"starbug2-match badinput.fits --output={TEST_PATH}") ==
            ExitStates.EXIT_EARLY)
    assert (run(f"starbug2-match badinput.txt --output={TEST_PATH}") ==
            ExitStates.EXIT_EARLY)
    starbug_main(
        f"starbug2 -D --output={TEST_PATH} {TEST_IMAGE_FITS} "
        f"{TEST_FILTER_STRING}".split())
    assert (run(f"starbug2-match --output={TEST_PATH} {IMAGE_AP_FITS}") ==
            ExitStates.EXIT_EARLY)


def test_match_basic_run_through():
    starbug_main(
        f"starbug2 --output={TEST_PATH} -Do {OUT_1_FITS}"
        f"   {TEST_IMAGE_FITS} -sPSF_FILE={TEST_PSF_FITS}"
        f" {TEST_FILTER_STRING}".split())
    starbug_main(
        f"starbug2 --output={TEST_PATH} -Do {OUT_2_FITS} "
        f" {TEST_IMAGE_FITS} -sPSF_FILE={TEST_PSF_FITS} "
        f" {TEST_FILTER_STRING}".split())
    assert (run(
        f"starbug2-match --output={TEST_PATH} -sPSF_FILE={TEST_PSF_FITS}"
        f" {OUT_1_AP_FITS}"
        f" {OUT_2_AP_FITS}"
        f" {TEST_FILTER_STRING}") == ExitStates.EXIT_SUCCESS)
    assert (run(
        f"starbug2-match"
        f" {OUT_1_AP_FITS}"
        f" {OUT_2_AP_FITS} -sPSF_FILE={TEST_PSF_FITS}"
        f" {TEST_FILTER_STRING} --output={TEST_PATH}") ==
            ExitStates.EXIT_SUCCESS)
    assert (run(
        f"starbug2-match"
        f" {OUT_1_AP_FITS}"
        f" {OUT_2_AP_FITS} -sPSF_FILE={TEST_PSF_FITS}"
        f" {TEST_FILTER_STRING} --output={TEST_PATH}") ==
            ExitStates.EXIT_SUCCESS)
    assert (run(
        f"starbug2-match"
        f" {OUT_1_AP_FITS}"
        f" {OUT_2_AP_FITS}"
        f" {TEST_FILTER_STRING} --output={TEST_PATH} "
        f"-sPSF_FILE={TEST_PSF_FITS}") ==
            ExitStates.EXIT_SUCCESS)
    assert (run(
        f"starbug2-match"
        f" {OUT_1_AP_FITS}"
        f" {OUT_2_AP_FITS}"
        f" {TEST_FILTER_STRING} --output={TEST_PATH} "
        f"-sPSF_FILE={TEST_PSF_FITS}") == ExitStates.EXIT_SUCCESS)
    assert (run(
        f"starbug2-match"
        f" {OUT_1_AP_FITS}"
        f" {OUT_2_AP_FITS}"
        f" {TEST_FILTER_STRING} -sPSF_FILE={TEST_PSF_FITS} "
        f"--output={TEST_PATH}") == ExitStates.EXIT_SUCCESS)


def test_mask():
    starbug_main(
        f"starbug2 --output={TEST_PATH} -Do "
        f"{OUT_1_FITS}  {TEST_IMAGE_FITS} -sPSF_FILE={TEST_PSF_FITS}"
        f" -s FILTER=F444W".split())
    starbug_main(
        f"starbug2 --output={TEST_PATH} -Do"
        f"{OUT_2_FITS}  {TEST_IMAGE_FITS} -sPSF_FILE={TEST_PSF_FITS}"
        f" -s FILTER=F444W".split())
    assert run(
        f"starbug2-match --output={TEST_PATH} -vmF444W>20 "
        f"{OUT_1_AP_FITS} {OUT_2_AP_FITS} "
        f"-sPSF_FILE={TEST_PSF_FITS}") == ExitStates.EXIT_SUCCESS


@pytest.fixture(autouse=True)
def init():
    clean()
