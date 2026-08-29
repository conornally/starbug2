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
from starbug2.constants import ExitStates
from tests.generic import (
    clean, TEST_IMAGE_FITS, TEST_FILTER_STRING, TEST_PATH_STR, TEST_PATH,
    TEST_FILTER_STRING_NO_G)

# different fit files paths
TEST_IMAGE_AP_FITS: Final[str] = os.path.join(
    TEST_PATH_STR, "image-ap.fits")
TEST_PSF_FITS: Final[str] = os.path.join(
    TEST_PATH_STR, "psf.fits")
TEST_IMAGE_BGD_FITS: Final[str] = os.path.join(
    TEST_PATH_STR, "image-bgd.fits")
TEST_IMAGE_RES_FIT: Final[str] = os.path.join(
    TEST_PATH_STR, "image-res.fits")
TEST_IMAGE_2_FITS: Final[str] = os.path.join(
    TEST_PATH_STR, "image2.fits")


def run(s):
    return starbug_main((s + f"  --output={TEST_PATH}").split())


def test_start():
    clean()
    assert run("starbug2 -h") == ExitStates.EXIT_EARLY
    assert run("starbug2 -vh") == ExitStates.EXIT_EARLY
    assert run("starbug2 --version") == ExitStates.EXIT_SUCCESS
    assert run("starbug2 -vDABPh") == ExitStates.EXIT_EARLY
    assert run("starbug2") == ExitStates.EXIT_FAIL
    clean()


def test_param():
    clean()
    assert run("starbug2 --local-param") == ExitStates.EXIT_SUCCESS
    assert run("starbug2 --update-param") == ExitStates.EXIT_SUCCESS
    assert (run(
        f"starbug2 -p starbug.param {TEST_IMAGE_FITS}"
        f" {TEST_FILTER_STRING_NO_G} "
        f"-sPSF_FILE={TEST_PSF_FITS} ") == ExitStates.EXIT_SUCCESS)
    clean()


def test_detect():
    clean()
    assert (run(
        f"starbug2 -v {TEST_IMAGE_FITS}"
        f" {TEST_FILTER_STRING_NO_G} -sPSF_FILE={TEST_PSF_FITS} ") ==
            ExitStates.EXIT_SUCCESS)
    assert (run(
        f"starbug2 -D {TEST_IMAGE_FITS} "
        f"{TEST_FILTER_STRING} -sPSF_FILE={TEST_PSF_FITS} ") ==
            ExitStates.EXIT_SUCCESS)
    assert (run(
        f"starbug2 --detect {TEST_IMAGE_FITS}"
        f" {TEST_FILTER_STRING} -sPSF_FILE={TEST_PSF_FITS} ") ==
            ExitStates.EXIT_SUCCESS)
    assert (run(f"starbug2 -D -sSIGSKY=3 -sSIGSRC=15 {TEST_IMAGE_FITS}"
                f" -sPSF_FILE={TEST_PSF_FITS}  {TEST_FILTER_STRING}") ==
            ExitStates.EXIT_SUCCESS)
    clean()


def test_bgd():
    clean()
    assert (run(
        f"starbug2 -D {TEST_IMAGE_FITS} "
        f"{TEST_FILTER_STRING} -sPSF_FILE={TEST_PSF_FITS} ") ==
            ExitStates.EXIT_SUCCESS)
    assert (run(f"starbug2 -d {TEST_IMAGE_AP_FITS}"
                f" -B {TEST_IMAGE_FITS} "
                f"{TEST_FILTER_STRING} -sPSF_FILE={TEST_PSF_FITS}") ==
            ExitStates.EXIT_SUCCESS)
    assert (run(f"starbug2 -d {TEST_IMAGE_AP_FITS} "
                f"--background {TEST_IMAGE_FITS} "
                f"{TEST_FILTER_STRING} -sPSF_FILE={TEST_PSF_FITS} ") ==
            ExitStates.EXIT_SUCCESS)
    assert (run(f"starbug2 -vf -B {TEST_IMAGE_FITS}"
                f" {TEST_FILTER_STRING} -sPSF_FILE={TEST_PSF_FITS} ") ==
            ExitStates.EXIT_SUCCESS)
    clean()


def test_psf():
    clean()
    assert (run(f"starbug2 -DB {TEST_IMAGE_FITS}"
                f" {TEST_FILTER_STRING} -sPSF_FILE={TEST_PSF_FITS} ") ==
            ExitStates.EXIT_SUCCESS)
    assert run(f"starbug2 -d {TEST_IMAGE_AP_FITS} -b "
               f"{TEST_IMAGE_BGD_FITS} -P {TEST_IMAGE_FITS} "
               f"-sPSF_FILE={TEST_PSF_FITS}"
               f" {TEST_FILTER_STRING}") == ExitStates.EXIT_SUCCESS
    assert run(f"starbug2 -fP {TEST_IMAGE_FITS} "
               f"-sPSF_FILE={TEST_PSF_FITS}"
               f" {TEST_FILTER_STRING}") == ExitStates.EXIT_SUCCESS
    assert run(f"starbug2 -d {TEST_IMAGE_AP_FITS}"
               f" -P {TEST_IMAGE_FITS} "
               f"-sPSF_FILE={TEST_PSF_FITS}"
               f" {TEST_FILTER_STRING}") == ExitStates.EXIT_SUCCESS
    assert run(f"starbug2 -fBP {TEST_IMAGE_FITS} "
               f"-sPSF_FILE={TEST_PSF_FITS}"
               f" {TEST_FILTER_STRING}") == ExitStates.EXIT_SUCCESS
    assert run(f"starbug2 -fPs GEN_RESIDUAL=1 {TEST_IMAGE_FITS} "
               f"-sPSF_FILE={TEST_PSF_FITS}"
               f" {TEST_FILTER_STRING}") == ExitStates.EXIT_SUCCESS
    clean()


def test_residual():
    clean()
    assert (run(f"starbug2 -DB {TEST_IMAGE_FITS} "
                f"{TEST_FILTER_STRING} -sPSF_FILE={TEST_PSF_FITS} ") ==
            ExitStates.EXIT_SUCCESS)
    assert run(
        f"starbug2 -sPSF_FILE={TEST_PSF_FITS} "
        f"-fSs GEN_RESIDUAL=1 {TEST_IMAGE_FITS}"
        f" {TEST_FILTER_STRING}") == ExitStates.EXIT_SUCCESS
    assert run(f"starbug2 {TEST_IMAGE_RES_FIT} "
               f"-sPSF_FILE={TEST_PSF_FITS}"
               f" {TEST_FILTER_STRING_NO_G}") == ExitStates.EXIT_SUCCESS
    assert run(f"starbug2 -D {TEST_IMAGE_RES_FIT} "
               f"-sPSF_FILE={TEST_PSF_FITS} "
               f" {TEST_FILTER_STRING}") == ExitStates.EXIT_SUCCESS
    assert run(f"starbug2 -fB {TEST_IMAGE_RES_FIT} "
               f"-sPSF_FILE={TEST_PSF_FITS} "
               f" {TEST_FILTER_STRING}") == ExitStates.EXIT_SUCCESS
    assert run(f"starbug2 -fP {TEST_IMAGE_RES_FIT} "
               f"-sPSF_FILE={TEST_PSF_FITS}"
               f" {TEST_FILTER_STRING}") == ExitStates.EXIT_SUCCESS
    assert run(f"starbug2 -fPs GEN_RESIDUAL=1 {TEST_IMAGE_RES_FIT} "
               f" -sPSF_FILE={TEST_PSF_FITS}"
               f" {TEST_FILTER_STRING}") == ExitStates.EXIT_SUCCESS

    assert run(f"starbug2 -fSA {TEST_IMAGE_FITS} "
               f"-sPSF_FILE={TEST_PSF_FITS} "
               f" {TEST_FILTER_STRING}") == ExitStates.EXIT_SUCCESS
    clean()


def test_n_cores():
    clean()
    os.system(f"cp {TEST_IMAGE_FITS} {TEST_IMAGE_2_FITS}")
    assert (run(
        f"starbug2 {TEST_IMAGE_FITS} {TEST_IMAGE_2_FITS}"
        f" {TEST_FILTER_STRING_NO_G} -sPSF_FILE={TEST_PSF_FITS} ") ==
            ExitStates.EXIT_SUCCESS)
    assert (run(
        f"starbug2 -n2 {TEST_IMAGE_FITS} -sPSF_FILE={TEST_PSF_FITS} "
        f"{TEST_IMAGE_2_FITS} {TEST_FILTER_STRING_NO_G}") ==
            ExitStates.EXIT_SUCCESS)
    assert (run(f"starbug2 -vD {TEST_IMAGE_FITS} {TEST_IMAGE_2_FITS}"
                f" {TEST_FILTER_STRING} -sPSF_FILE={TEST_PSF_FITS} ") ==
            ExitStates.EXIT_SUCCESS)

    with pytest.raises(
            ValueError,
            match="Number of processes must be at least 1"):
        run(f"starbug2 -Dn0 {TEST_IMAGE_FITS} {TEST_IMAGE_2_FITS}"
            f" {TEST_FILTER_STRING} -sPSF_FILE={TEST_PSF_FITS}")
    assert (run(f"starbug2 -Dn1 {TEST_IMAGE_FITS} {TEST_IMAGE_2_FITS}"
                f" {TEST_FILTER_STRING} -sPSF_FILE={TEST_PSF_FITS} ") ==
            ExitStates.EXIT_SUCCESS)
    assert (run(f"starbug2 -Dn2 {TEST_IMAGE_FITS} {TEST_IMAGE_2_FITS}"
                f" {TEST_FILTER_STRING} -sPSF_FILE={TEST_PSF_FITS} ") ==
            ExitStates.EXIT_SUCCESS)
    assert (run(f"starbug2 -Dn4 {TEST_IMAGE_FITS} {TEST_IMAGE_2_FITS}"
                f" {TEST_FILTER_STRING} -sPSF_FILE={TEST_PSF_FITS} ") ==
            ExitStates.EXIT_SUCCESS)

    assert run(f"starbug2 -DBP {TEST_IMAGE_FITS} {TEST_IMAGE_2_FITS} "
               f"-sPSF_FILE={TEST_PSF_FITS} "
               f" {TEST_FILTER_STRING}") == ExitStates.EXIT_SUCCESS
    assert run(f"starbug2 -vDBPn2 {TEST_IMAGE_FITS} {TEST_IMAGE_2_FITS}"
               f" -sPSF_FILE={TEST_PSF_FITS}"
               f" {TEST_FILTER_STRING}") == ExitStates.EXIT_SUCCESS
    assert (run(f"starbug2 -DM {TEST_IMAGE_FITS} {TEST_IMAGE_2_FITS}"
                f" {TEST_FILTER_STRING} -sPSF_FILE={TEST_PSF_FITS} ") ==
            ExitStates.EXIT_SUCCESS)
    assert (run(f"starbug2 -DMn2 {TEST_IMAGE_FITS} "
                f"{TEST_IMAGE_2_FITS} {TEST_FILTER_STRING} "
                f"-sPSF_FILE={TEST_PSF_FITS} ") ==
            ExitStates.EXIT_SUCCESS)

    assert (run(f"starbug2 -D {TEST_IMAGE_AP_FITS} "
                f"{TEST_IMAGE_FITS} {TEST_FILTER_STRING} "
                f"-sPSF_FILE={TEST_PSF_FITS} ") ==
            ExitStates.EXIT_MIXED)
    assert run(f"starbug2 -D bad.fits {TEST_IMAGE_FITS}"
               f" {TEST_FILTER_STRING} "
               f"-sPSF_FILE={TEST_PSF_FITS} ") == ExitStates.EXIT_MIXED
    clean()
