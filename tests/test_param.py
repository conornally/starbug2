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
from starbug2.constants import STAR_BUG_PARAMS, PROBLEMATIC_FILTER_WARNING
from starbug2.core.star_bug_config import StarBugMainConfig
from tests.generic import TEST_PATH_STR, verify_test_data_exists, TEST_PATH

TEST_PARAM_PATH: Final[str] = os.path.join(
    TEST_PATH_STR, "../param_files/old_format.param")


def test_parse_param():
    verify_test_data_exists()
    assert type(StarBugMainConfig.parse_param("A=1//.")) == dict

    assert StarBugMainConfig.parse_param("A = 1 //.") == {'A': 1}
    assert StarBugMainConfig.parse_param("A = B //.") == {'A': 'B'}
    assert StarBugMainConfig.parse_param("A = B //.\n") == {'A': 'B'}

    assert StarBugMainConfig.parse_param("A = //.\n") == {'A': ''}
    assert StarBugMainConfig.parse_param(" = //.") == {}

    assert StarBugMainConfig.parse_param("A=B") == {"A": "B"}
    assert StarBugMainConfig.parse_param("A=B/") == {"A": "B/"}
    assert StarBugMainConfig.parse_param("A=B/.") == {"A": "B/."}
    assert StarBugMainConfig.parse_param("A=1/.") == {"A": "1/."}

    assert StarBugMainConfig.parse_param("A      =1") == {"A": 1}
    assert StarBugMainConfig.parse_param("A=1      ") == {"A": 1}
    assert StarBugMainConfig.parse_param("A=1     a") == {"A": "1     a"}


def test_load_default_params():
    config: StarBugMainConfig = StarBugMainConfig()
    assert config.param_tag == STAR_BUG_PARAMS


def test_load_params():
    config: StarBugMainConfig = StarBugMainConfig.load_params("does_not_exist")

    os.system(f"starbug2 --local-param --output={TEST_PATH}")
    second_config: StarBugMainConfig = (
        StarBugMainConfig.load_params("starbug.param"))

    for value, _ in config.MAIN_PARAM_FILE_MAP.values():
        assert getattr(config, value) == getattr(second_config, value)
    assert second_config.param_tag == STAR_BUG_PARAMS
    assert config.param_tag == STAR_BUG_PARAMS
    os.remove(os.path.join(TEST_PATH, "starbug.param"))


def test_update_params():
    os.system(f"starbug2 --local-param --output={TEST_PATH}")
    os.system(
        f"sed -i s/PARAM/PARAM1/g {os.path.join(TEST_PATH, 'starbug.param')}")

    with pytest.raises(
        TypeError,
        match="Param PARAM1 no longer works within Starbug2. Please "
              "execute starbug2 --update-param"):
        StarBugMainConfig.load_params(os.path.join(TEST_PATH, 'starbug.param'))
    os.remove(os.path.join(TEST_PATH, 'starbug.param'))


def test_update_params_no_changes():
    os.system(f"starbug2 --local-param --output={TEST_PATH}")
    os.system(f"starbug2 --update-param --output={TEST_PATH}")
    os.remove(os.path.join(TEST_PATH, 'starbug.param'))


def test_update_params_old_to_new():
    os.system(f"starbug2 -p {TEST_PARAM_PATH} --update-param")


def test_f150w2_filter(capsys):
    config: StarBugMainConfig = StarBugMainConfig()
    config.custom_filter = "F150W2"

    # check for warning
    captured = capsys.readouterr()
    assert PROBLEMATIC_FILTER_WARNING in captured.err
