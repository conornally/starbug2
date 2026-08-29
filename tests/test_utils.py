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

import numpy as np
from astropy.table import Table
from astropy.io import fits

from starbug2.constants import Units
from starbug2.utilities import utils
from tests.generic import check_shape


def test_str_nk_tn() -> None:
    assert utils.append_chars("", 3, 'a') == "aaa"
    assert utils.append_chars("a", 3, 'a') == "aaaa"
    assert utils.append_chars("a", 0, 'a') == "a"


def test_split_f_name() -> None:
    f_name = "/path/to/file.fits"
    d, f, e = utils.split_file_name(f_name)
    assert d == "/path/to"
    assert f == "file"
    assert e == ".fits"

    f_name = "file.fits"
    d, f, e = utils.split_file_name(f_name)
    assert d == "."
    assert f == "file"
    assert e == ".fits"

    f_name = "file"
    d, f, e = utils.split_file_name(f_name)
    assert d == "."
    assert f == "file"
    assert e == ""


def test_flux2mag() -> None:
    # Input shape validation
    assert len(utils.flux_to_pogson_mag(1, None)[0]) == 1
    assert len(utils.flux_to_pogson_mag(np.ones(10), None)[0]) == 10
    assert len(utils.flux_to_pogson_mag(np.full(10, np.nan), None)[0]) == 10

    a, b = utils.flux_to_pogson_mag(np.empty(10), np.empty(10))
    assert len(a) == len(b)
    a, b = utils.flux_to_pogson_mag(1, 1)
    assert len(a) == len(b)
    a, b = utils.flux_to_pogson_mag(0, 0)
    assert len(a) == len(b)

    # Normal flux validation
    flux = np.array([1, 100, 999, 123, 3.4, 87654, np.pi])
    flux_err = None
    mag, mag_err = utils.flux_to_pogson_mag(flux, flux_err)
    assert np.all(np.isclose(mag, -2.5 * np.log10(flux)))

    # Boundary fluxes
    flux = np.array([0, 0.0, -1, np.nan])
    flux_err = None
    mag, mag_err = utils.flux_to_pogson_mag(flux, flux_err)
    assert np.isnan(mag).all()

    # should be -inf
    assert utils.flux_to_pogson_mag(np.inf)[0] == -np.inf

    # Should be nan
    assert np.isnan(utils.flux_to_pogson_mag(-np.inf)[0])

    # flux_err
    flux = np.array([1234, 1, 0.00001, 10])
    flux_err = np.array([1, 100, 123456, 1.234567])
    mag, mag_err = utils.flux_to_pogson_mag(np.ones(flux.shape), flux_err)

    # flux all 1
    assert np.all(np.equal(
        mag_err, 2.5 * np.log10(1.0 + (flux_err / np.ones(flux.shape)))))
    mag, mag_err = utils.flux_to_pogson_mag(flux, flux_err)

    # random fluxes
    assert np.all(np.equal(
        mag_err, 2.5 * np.log10(1.0 + (flux_err / flux))))

    # boundary flux_errs
    assert utils.flux_to_pogson_mag(1, None)[1] == 0
    assert np.isnan(utils.flux_to_pogson_mag(1, np.nan)[1])
    assert np.isnan(utils.flux_to_pogson_mag(1, -1)[1])


def test_find_col_names() -> None:
    # noinspection SpellCheckingInspection
    tab = Table(
        None, names=["A", "word", "word1", "word2", "notword", "_word"])
    res = utils.find_col_names(tab, "word")

    assert res is not None
    assert res == ["word", "word1", "word2"]
    # noinspection SpellCheckingInspection
    assert utils.find_col_names(tab, "badmatch") == []


def test_tab_append() -> None:
    base = Table([[0, 0], [0, 0]], names=('a', 'b'))
    tab = Table([[1, 1], [1, 1]], names=('a', 'b'))
    exp = Table([[0, 0, 1, 1], [0, 0, 1, 1]], names=('a', 'b'))
    out = utils.combine_tables(base, tab)

    # Safely compare tables via their underlying numpy structures
    assert np.all(out.as_array() == exp.as_array())

    tab1 = tab.copy()
    out = utils.combine_tables(None, tab1)
    assert np.all(out.as_array() == tab.as_array())


def test_parse_unit() -> None:
    assert utils.parse_unit("10p") == (10, Units.PIX)
    assert utils.parse_unit("10s") == (10, Units.ARCSEC)
    assert utils.parse_unit("10m") == (10, Units.ARCMIN)
    assert utils.parse_unit("10d") == (10, Units.DEG)

    assert utils.parse_unit("10.1s") == (10.1, Units.ARCSEC)
    assert utils.parse_unit("-10.1s") == (-10.1, Units.ARCSEC)
    assert utils.parse_unit("0s") == (0, Units.ARCSEC)
    assert utils.parse_unit("0") == (0, None)

    assert utils.parse_unit("") == (None, None)
    assert utils.parse_unit("p") == (None, None)


def test_remove_duplicates() -> None:
    lst = ["a", "b", "b", "c", "b", "c"]
    lst2 = utils.remove_duplicates(lst)
    assert lst2 == ["a", "b", "c"]

    assert utils.remove_duplicates([]) == []
    assert utils.remove_duplicates(["a"]) == ["a"]


def test_h_cascade() -> None:
    t1 = [[1, 1, 0],
          [2, 2, 0],
          [3, 3, 0]]

    t2 = [[1, 1, 0],
          [2, 2, 0],
          [3, 3, 1],
          [4, 4, 0]]

    tables = [
        Table(np.array(t1), names=["A", "B", "flag"],
              dtype=[float, float, np.uint16]),
        Table(np.array(t2), names=["A", "B", "flag"],
              dtype=[float, float, np.uint16])
    ]
    nan = np.nan
    res = utils.h_cascade(tables)

    # Corrected alignment: Since t1 has fewer rows than t2, missing slots
    # belong to t1 components
    test = Table(
        np.ma.array([
            [1, 1, 0, 1, 1, 0],
            [2, 2, 0, 2, 2, 0],
            [3, 3, 0, 3, 3, 1],
            [4, 4, 0, nan, nan, 0]
        ]),
        dtype=[float, float, np.uint16, float, float, np.uint16],
        names=["A_1", "B_1", "flag_1", "A_2", "B_2", "flag_2"]
    )

    res = utils.fill_nan(res)
    check_shape(res, test)


def test_collapse_header() -> None:
    # noinspection SpellCheckingInspection
    header = fits.Header({
        "OK": 0,
        "PARAMFILE": "/PATH/TO/FILE/THAT/IS/TOO/LONG/FOR/A/HIERARCH/CARD",
        "PARAMFILE2": "/PATH/TO/FILE/THAT/IS/TOO/LONG/FOR/A/HIERARCH/CARD"
    })

    h = utils.collapse_header(header)
    assert h["COMMENT"] is not None
    # Explicit type validation using isinstance instead of un-idiomatic direct
    # type comparison
    assert isinstance(utils.collapse_header({"a": "b"}), fits.Header)
