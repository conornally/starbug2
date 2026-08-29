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
from typing import override, Final, Any

import numpy as np
import astropy.units as u
from astropy.table import Table, hstack
from starbug2.constants import HeaderTags, TableColumn
from starbug2.utilities.filters import STAR_BUG_FILTERS
from starbug2.matching.generic_match import GenericMatch
from starbug2.utilities.utils import (
    Loading, printf, p_error, fill_nan, find_col_names, warn, puts)

# keys for catalogue fields.
# noinspection SpellCheckingInspection
_OBS: Final[str] = "OBSERVTN"
_VISIT: Final[str] = "VISIT"
_EXPOSURE: Final[str] = "EXPOSURE"
_DEFAULT: Final[str] = "DEFAULT"

# Hard coding separations for now
_SEPARATION_VALUES: Final[dict[str, float]] = {
    "F277W": 0.1,
    "F560W": 0.15,
    "F1000W": 0.2,
    "F1500W": 0.25,
    _DEFAULT: 0.06,
}


class BandMatch(GenericMatch):
    # filter flag for kwargs.
    # noinspection SpellCheckingInspection
    FILTER: Final[str] = "fltr"
    THRESHOLD: Final[str] = "threshold"

    # match methods
    _FIRST: Final[str] = "first"
    _LAST: Final[str] = "last"
    _BOOT_STRAP: Final[str] = "bootstrap"

    # warning messages
    _WRONG_THRESHOLD: Final[str] = (
        "Threshold values must be scalar or list with length 1 less than the "
        "catalogue list. The final element is being ignored.\n")

    def __init__(self, **kwargs: Any) -> None:
        if BandMatch.FILTER in kwargs:
            if not isinstance(kwargs[BandMatch.FILTER], list):
                warn(f"{self.FILTER} input should be a list, "
                     "there may be unexpected behaviour\n")
                self._filter_list: list[str] = [kwargs[BandMatch.FILTER]]
            else:
                self._filter_list: list[str] = kwargs[BandMatch.FILTER]
            del kwargs[BandMatch.FILTER]
        else:
            self._filter_list: list[str] = []

        if BandMatch.THRESHOLD in kwargs:
            if isinstance(kwargs[BandMatch.THRESHOLD], list):
                kwargs[BandMatch.THRESHOLD] = (
                    np.array(kwargs[BandMatch.THRESHOLD], dtype=object))

        super().__init__(**kwargs, method="Band Matching")

    @property
    def filter_list(self) -> list[str]:
        return self._filter_list

    def order_catalogues(self, catalogues: list[Table]) -> list[Table]:
        """
        Reorder catalogue list into increasing wavelength size.
        This only works for JWST bands. Unrecognised filters will be left
        unchanged. The function should also set the self.FILTER variable if
        possible.

        :param catalogues: List of astropy tables with meta keys 'FILTER'
        :type catalogues: list[astropy.table.Table]
        :return: The same list reordered by wavelength
        :rtype: list[astropy.table.Table]
        """

        status: int = -1
        _ii = None
        filter_list: list[str] = self.filter_list
        if filter_list is None:
            raise Exception("no filer was set")

        sorters = [
            # META in JWST filters
            lambda t: list(STAR_BUG_FILTERS.keys()).index(
                t.meta.get(HeaderTags.FILTER)),

            # col_names in JWST filters
            lambda t: list(STAR_BUG_FILTERS.keys()).index(
                (set(t.colnames) & set(STAR_BUG_FILTERS.keys())).pop()),

            # META in self.filters
            lambda t: filter_list.index(t.meta.get(HeaderTags.FILTER)),

            # col_names in JWST filters
            lambda t: filter_list.index(
                (set(t.colnames) & set(filter_list)).pop())
        ]

        for n, fn in enumerate(sorters):
            try:
                catalogues.sort(key=fn)
                _ii = map(fn, catalogues)
                status = n
                break
            except (KeyError, AttributeError, TypeError, ValueError) as e:
                warn(f"failed to use sorter {n} due to {str(e)}\n")
                pass

        if status < 0:
            p_error(
                "Unable to reorder catalogues, leaving input order"
                " untouched.\n")
        elif status <= 1 and (_ii is not None):
            # JWST filters
            self._filter_list = [list(STAR_BUG_FILTERS.keys())[i] for i in _ii]

        self._load: Loading = Loading(sum(len(c) for c in catalogues[1:]))

        return catalogues

    def jwst_order(self, catalogues: list[Table]):
        pass

    @override
    def match(
            self, catalogues: list[Table], method: str = "first",
            **kwargs: Any) -> Table:
        # noinspection SpellCheckingInspection
        """
        Given a list of catalogues, it will reorder them into increasing
        wavelength or to match the fltr= keyword in the initialiser.
        The matching then uses the shortest wavelength available position.
        I.e. If F115W, F444W, F770W are input, the F115W centroid positions
        will be taken as "correct". If a source is not resolved in this band,
        the next most astrometric ally accurate position is taken, i.e. F444W

        :param catalogues: List of `astropy.table.Table` objects containing
        the meta item "FILTER=XXX"
        :type catalogues: List of `astropy.table.Table`
        :param method: Centroid method
            "first" -   Use the position corresponding to the earliest
                        appearance of the source
            "last"  -   Use the position corresponding to the latest
                        appearance of the source
            "bootsrap"- ???
            "average" - ???
        :type method: str
        :param kwargs:
        :return: Matched catalogue
        :rtype: astropy.Table
        """
        catalogues = self.order_catalogues(catalogues)

        printf("Bands: %s\n" % ', '.join(self.filter_list))

        assert self.threshold is not None
        assert len(self.filter_list) != 0
        if type(self.threshold) in (list, np.ndarray):
            if len(self.threshold) != (len(catalogues) - 1):
                warn(self._WRONG_THRESHOLD)
                self._threshold = self._threshold[:-1]
        else:
            self._threshold = (
                np.full(len(catalogues) - 1, self._threshold) * u.arcsec)

        assert self._threshold is not None
        threshold_strs: list[str] = []
        for threshold in self.threshold:
            if hasattr(threshold, 'value'):
                threshold_strs.append(f"{threshold.value}")
            else:
                threshold_strs.append(f"{threshold}")
        printf(f"Thresholds: {', '.join(threshold_strs)}\n")

        if self._col_names is None:
            self._col_names = [
                TableColumn.RA, TableColumn.DEC, TableColumn.FLAG,
                TableColumn.NUM,
                *self._filter_list, *["e%s" % f for f in self._filter_list]]

        assert self._col_names is not None
        printf("Columns: %s\n" % ", ".join(self._col_names))

        if method not in (self._FIRST, self._LAST, self._BOOT_STRAP):
            method = self._FIRST

        # Begin
        base: Table = self.build_meta(catalogues)
        _threshold = self._threshold.copy()
        for n, tab in enumerate(catalogues):
            # Temporarily recast threshold
            self._threshold = _threshold[n - 1]
            assert self._threshold is not None
            self._load.msg = (
                f"{self.filter_list[n]} ("
                f"{np.array2string(self._threshold)}\")")
            col_names = [
                name for name in self._col_names if name in tab.colnames]

            tmp = self.inner_match(base, tab, join_type="or")

            col_names.remove(TableColumn.RA)
            col_names.remove(TableColumn.DEC)
            base = fill_nan(hstack((base, tmp[col_names])))
            base.rename_columns(
                col_names, ["%s_%d" % (name, n + 1) for name in col_names])

            if TableColumn.RA not in base.colnames:
                base = fill_nan(hstack(
                    (tmp[TableColumn.RA, TableColumn.DEC], base)))
            elif method == self._FIRST:
                _mask = np.logical_and(np.isnan(
                    base[TableColumn.RA]), tmp[TableColumn.RA] != np.nan)
                base[TableColumn.RA][_mask] = tmp[TableColumn.RA][_mask]
                base[TableColumn.DEC][_mask] = tmp[TableColumn.DEC][_mask]
            elif method == self._LAST:
                _mask = ~np.isnan(tmp[TableColumn.RA])
                base[TableColumn.RA][_mask] = tmp[TableColumn.RA][_mask]
                base[TableColumn.DEC][_mask] = tmp[TableColumn.DEC][_mask]
            elif method == self._BOOT_STRAP:
                _mask = ~np.isnan(tmp[TableColumn.RA])
                base.rename_columns(
                    (TableColumn.RA, TableColumn.DEC),
                    ("_RA_%d" % n, "_DEC_%d" % n))
                base = hstack((base, tmp[[TableColumn.RA, TableColumn.DEC]]))

        # Set threshold back at the end
        self._threshold = _threshold

        # Fix column names
        for name in self._col_names:
            all_cols = find_col_names(base, name)
            if len(all_cols) == 1:
                base.rename_column(all_cols.pop(), name)

        # Finalise NUM and flag column
        tmp = self.finish_matching(base, col_names=["NUM", "flag"])
        base.remove_columns(
            (*find_col_names(base, "NUM"), *find_col_names(base, "flag")))
        base.add_column(tmp["NUM"], index=2)
        base.add_column(tmp["flag"], index=3)

        return base

    def band_match(
            self, catalogues, col_names=(TableColumn.RA, TableColumn.DEC)):
        """
        Given a list of catalogues (with filter names in the metadata), match
        them in order of decreasing astrometric accuracy. If F115W, F444W,
        F770W are input, the F115W centroid positions will be taken as
        "correct". If a source is not resolved in this band, the next most
        astrometric ally accurate position is taken, i.e. F444W

        :param catalogues: list of tables
        :param col_names: the col names to match against.
        :return: Matched catalogue
        :rtype: astropy.Table
        """

        # ORDER the tables into the correct order (increasing wavelength)
        tables = np.full(len(STAR_BUG_FILTERS), None)
        mask = np.full(len(STAR_BUG_FILTERS), False)
        for tab in catalogues:
            if HeaderTags.FILTER in tab.meta.keys():
                if tab.meta[HeaderTags.FILTER] in STAR_BUG_FILTERS:
                    ii = list(STAR_BUG_FILTERS.keys()).index(
                        tab.meta[HeaderTags.FILTER])
                    tables[ii] = tab
                    mask[ii] = True
                else:
                    p_error(
                        "Unknown filter '%s' (skipping)..\n" %
                        tab.meta[HeaderTags.FILTER])
            elif _tmp := set(STAR_BUG_FILTERS.keys()) & set(tab.col_names):
                ii = list(STAR_BUG_FILTERS.keys()).index(_tmp.pop())
                tables[ii] = tab
                mask[ii] = True
            else:
                p_error("Cannot find 'FILTER' in table meta (skipping)..\n")

        # document bands
        s = "Bands: "
        for filter_string, tab in zip(STAR_BUG_FILTERS.keys(), tables):
            if tab:
                s += "%5s " % filter_string
        puts(s)

        # Match in increasing wavelength order
        base = Table(None)
        load = Loading(
            sum([len(t) for t in tables[mask][1:]]), "matching", res=100)
        for filter_string, tab in zip(STAR_BUG_FILTERS.keys(), tables):
            if not tab:
                continue

            # removing empty magnitude rows
            tab.remove_rows(np.isnan(tab[filter_string]))
            load.msg = "matching:%s" % filter_string
            _col_names = (
                list(name for name in tab.col_names if name in col_names))
            if not len(base):
                tmp = tab[_col_names].copy()
            else:
                idx, d2d, _ = self.inner_match(base, tab)
                tmp = Table(
                    np.full((len(base), len(_col_names)), np.nan),
                    names=_col_names)

                # Hard coding separations for now
                separation = _SEPARATION_VALUES.get(
                    filter_string, _SEPARATION_VALUES[_DEFAULT])

                for ii, (src, IDX, sep) in enumerate(zip(tab, idx, d2d)):
                    load.msg = (
                        "matching:%s(%.2g\")" % (filter_string, separation))
                    load()
                    load.show()

                    if ((sep <= separation * u.arcsec)
                            and (sep == min(d2d[idx == IDX]))):
                        for name in _col_names:
                            tmp[IDX][name] = src[name]
                    else:
                        tmp.add_row(src[_col_names])

            tmp.rename_column(TableColumn.FLAG, "flag_%s" % filter_string)
            base = hstack((
                base, tmp[[filter_string, "e%s" % filter_string,
                           "flag_%s" % filter_string]]
            ))

            base = (Table(
                base, dtype=[float] * len(base.col_names))
                    .filled(np.nan))  # type: ignore

            # Only keep the most astronomically correct position
            if TableColumn.RA not in base.col_names:
                base = hstack((tmp[[TableColumn.RA, TableColumn.DEC]], base))
            else:
                _mask = np.logical_and(
                    np.isnan(base[TableColumn.RA]),
                    tmp[TableColumn.RA] != np.nan)
                base[TableColumn.RA][_mask] = tmp[TableColumn.RA][_mask]
                base[TableColumn.DEC][_mask] = tmp[TableColumn.DEC][_mask]

        # Sort out flags
        flag: np.ndarray = np.zeros(len(base), dtype=np.uint16)
        for f_col in find_col_names(base, TableColumn.FLAG):
            flag |= base[f_col].value.astype(np.uint16)
            base.remove_column(f_col)
        base.add_column(flag, name=TableColumn.FLAG)

        return base.filled(np.nan)  # type: ignore
