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

import glob
import os
import pytest

from starbug2.utilities.filters import STAR_BUG_FILTERS
from starbug2.constants import STARBUG_DATA_DIR
from starbug2.jwst_support.initialise_psf_data import init_starbug_for_jwst
from starbug2.core.star_bug_config import StarBugMainConfig
from tests.generic import verify_test_data_exists


@pytest.mark.skipif(
    os.getenv("RUN_STAR_BUG_PRODUCTION_TESTS") is None or
    os.getenv("RUN_STAR_BUG_PRODUCTION_TESTS") == "false",
    reason="test_init locked out of normal development runs due to "
           "length of time to run, CPU resources required which nearly slags"
           " the machine."
)
def test_init():
    verify_test_data_exists()
    os.environ[STARBUG_DATA_DIR] = "/tmp/starbug"
    d = os.getenv(STARBUG_DATA_DIR)
    init_starbug_for_jwst(StarBugMainConfig())

    for f in STAR_BUG_FILTERS:
        assert glob.glob("%s/*%s*" % (d, f))
