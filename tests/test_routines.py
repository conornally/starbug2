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

from astropy.io import fits
from astropy.table import Table

from starbug2.constants import SCI, TableColumn
from starbug2.routines.background_estimate_routine import (
    BackGroundEstimateRoutine)
from starbug2.routines.detection_routines import DetectionRoutine
from starbug2.routines.source_properties import SourceProperties
from tests.generic import TEST_IMAGE_FITS, verify_test_data_exists


class TestDetection:
    verify_test_data_exists()
    im = fits.open(TEST_IMAGE_FITS)[SCI].data
    a = Table([[0, 10], [0, 10]],
              names=[TableColumn.X_CENTROID, TableColumn.Y_CENTROID])
    b = Table([[20, 10, 50], [20, 10, 0]],
              names=[TableColumn.X_CENTROID, TableColumn.Y_CENTROID])

    def test_detection_routine_none(self):
        dt = DetectionRoutine()
        assert len(dt.find_stars(None)) == 0

    def test_detection_routine_crashes(self):
        dt = DetectionRoutine()
        out = dt.find_stars(self.im.copy())
        assert len(out) != 0

    def test_detection_match(self):
        dt = DetectionRoutine()
        _a = self.a.copy()
        _b = self.b.copy()
        c = dt.match(_a, _b)
        assert isinstance(c, Table)
        assert len(_a) == len(self.a)
        assert len(_b) == len(self.b)
        assert len(c) == 4

    def test_bkg2d(self):
        b = DetectionRoutine().bkg2d(self.im.copy())
        assert type(b) is type(self.im)
        assert b.shape == self.im.shape


class TestBackground:
    def test_background_estimate_routine_none(self):
        bg = BackGroundEstimateRoutine(None)
        assert bg(None) is None


def test_source_properties_none():
    sp = SourceProperties(None, None)
    assert sp.calculate_crowding() is None
    assert sp.calculate_geometry(1) is None
