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
from typing import Any
import numpy as np
from astropy import units
from astropy.units.quantity import Quantity
from astropy.coordinates import SkyCoord
from astropy.table import Table, hstack, Column, vstack
from starbug2.constants import HeaderTags, SourceFlags, TableColumn
from starbug2.core.star_bug_config import StarBugMainConfig
from starbug2.utilities.utils import (
    Loading, printf, remove_duplicates, p_error, fill_nan, tab2array,
    find_col_names, flux_to_pogson_mag)


class GenericMatch:

    @staticmethod
    def build_meta(catalogues: list[Table]) -> Table:
        """
        Extracts structural tracking headers to set up baseline combined
        outputs.
        """
        meta: dict[str, Any] = catalogues[0].meta
        base: Table = Table(None, meta=meta)
        return base

    @staticmethod
    def mask_catalogues(
            catalogues: list[Table],
            mask: list[np.ndarray |
                       list[Any] | None] | np.ndarray | None) -> Table:
        """ Takes catalogues and masks and removes catalogues which don't
        match the mask.

        :param catalogues: the catalogues to match.
        :type catalogues: list (astropy.Tables)
        :param mask: the mask to apply.
        :return: an astro table with masked catalogues.
        :rtype: astropy.Table
        """
        masked: Table = Table(None)

        if mask is None or not isinstance(mask, (list, np.ndarray)):
            return masked
        if len(catalogues) != len(mask):
            return masked

        for subset, cat in zip(mask, catalogues):
            if subset is not None:
                if isinstance(subset, list):
                    subset = np.array(subset)

                if len(subset) == len(cat):
                    # 1. Grab rows where mask is False (the ones we want to
                    # remove)
                    masked_rows = cat[~subset]
                    masked = vstack((masked, masked_rows))

                    # 2. Extract indices where mask is False to safely strip
                    # them out
                    indices_to_remove: np.ndarray = np.where(~subset)[0]
                    cat.remove_rows(indices_to_remove.tolist())
        return masked

    @staticmethod
    def _sky_coords_not_cartesian(base: Table) -> SkyCoord:
        """
        create sky coords which are not cartesian

        :param base: the base.
        :type base: astropy.table.Table
        :return: a sky coords sing base values.
        :rtype: SkyCoord.
        """
        ra_cols: list[str] = list(
            name for name in base.colnames if TableColumn.RA in name)
        dec_cols: list[str] = list(
            name for name in base.colnames if TableColumn.DEC in name)
        ra: np.ndarray = np.nanmean(
            tab2array(base, col_names=ra_cols), axis=1)
        dec: np.ndarray = np.nanmean(
            tab2array(base, col_names=dec_cols), axis=1)
        return SkyCoord(ra=ra * units.deg, dec=dec * units.deg)

    @staticmethod
    def _sky_coords_cartesian(base: Table) -> SkyCoord:
        """
        create sky coords which are cartesian

        :param base: the base.
        :type base: astropy.table.Table
        :return: a sky coords sing base values.
        :rtype: SkyCoord.
        """
        x_cols: list[str] = list(
            name for name in base.colnames if name[0] == "x")
        y_cols: list[str] = list(
            name for name in base.colnames if name[0] == "y")
        x: np.ndarray = np.nanmean(tab2array(base, col_names=x_cols), axis=1)
        y: np.ndarray = np.nanmean(tab2array(base, col_names=y_cols), axis=1)
        return SkyCoord(
            x=x, y=y, z=np.zeros(len(x)), representation_type="cartesian")

    def __init__(
        self,
        threshold: Quantity | np.ndarray | None = None,
        col_names: list[str] | None = None,
        filter_string: str | None = None,
        verbose: int | None = None,
        p_file: str | None = None,
        method: str = "Generic Matching",
        load: Loading = Loading(1)
    ) -> None:
        """
        constructor for the generic match.

        :param threshold: Separation threshold in arc-seconds
        :type threshold: float or None
        :param col_names: List of str column names to include in the matching.
        Everything else will be discarded.
        :type col_names: list or None
        :param filter_string: Specifically set the filter of the catalogues
        :type filter_string: str or None
        :param verbose: Include verbose outputs
        :type verbose: int or None
        :param p_file: Parameter filename
        :type p_file: str or None
        :param load: the loading object
        :type load: Loading
        """
        config = StarBugMainConfig.load_params(p_file)

        self._threshold: Quantity | np.ndarray | None = threshold
        self._filter: str | None = config.custom_filter
        self._verbose: int | None = config.verbose_logs
        self.method: str = method

        if filter_string is not None:
            self._filter = filter_string
        if verbose is not None:
            self._verbose = verbose

        self._col_names: list[str] | None = col_names
        self._load: Loading = load

    def log(self, msg: str) -> None:
        """
        logs messages only when in verbose mode.
        :param msg: message to log
        :return: None
        """
        if self._verbose:
            printf(msg)

    def __str__(self) -> str:
        """
        string representation fo the generic match class.
        :return: str
        """
        threshold_string: str
        if isinstance(self._threshold, np.ndarray):
            threshold_string = np.array2string(self._threshold)
        else:
            threshold_string = f"{self._threshold}"

        s: list[str] = [
            f"{self.method}:",
            f"Filter: {self._filter}",
            f"Col names: {self._col_names}",
            f'Threshold: {threshold_string}"'
        ]
        return "\n".join(s)

    def __call__(self, *args: Any, **kwargs: Any) -> Table:
        """
        main entrance method.

        :param args: args to the matching algorithm.
        :param kwargs: extra args.
        :return: matched catalogue
        :rtype: astropy.Table
        """
        return self.match(*args, **kwargs)

    def init_catalogues(self, catalogues: list[Table]) -> list[Table]:
        # noinspection SpellCheckingInspection
        """
        This function is a bit of a "do everything" function

        It takes the input catalogues and removes any columns that aren't
        included in the self.colnames list. If this is None, then all columns
        are kept. Additionally, it initialises the loading bar with the summed
        length of all the input catalogues.
        Finally, it attempts to set the photometric filter that is being used

        :param catalogues: The input catalogues to work on
        :type catalogues: list (astropy.Tables)
        :return: The cleaned list of input catalogues
        :rtype: list (astropy.Table)
        """
        if len(catalogues) >= 2:
            self._load = Loading(
                sum(len(cat) for cat in catalogues[1:]), msg="initialising")
            if self._verbose:
                self._load.show()

        # initialise the column names if it wasn't already set
        col_names: list[str]
        if self._col_names is None:
            col_names = []
            for cat in catalogues:
                col_names += cat.colnames
            self._col_names = col_names
        else:
            col_names = self._col_names

        self._col_names = remove_duplicates(col_names)

        # clean out the column names not included in self._col_names
        for n, catalogue in enumerate(catalogues):
            keep: list[str] = list(
                set(catalogue.colnames) & set(col_names))
            keep = sorted(
                keep,
                key=lambda s:
                    self._col_names.index(s) if self._col_names else 0)
            catalogues[n] = catalogue[keep]

        if not self._filter:
            filter_string: str | None = (
                catalogues[0].meta.get(HeaderTags.FILTER))
            if filter_string is None:
                filter_string = "MAG"
            self._filter = filter_string

        return catalogues

    def match(
            self,
            catalogues: list[Table],
            join_type: str = "or",
            mask: (
                list[np.ndarray | list[Any] | None]
                | np.ndarray | None) = None,
            cartesian: bool = False,
            **kwargs: Any) -> Table:
        """
        This matching works as a basic match. Everything is included and the
        column names have _N appended to the end.

        :param catalogues: the catalogues to match
        :type catalogues: list (astropy.Tables)
        :param join_type: Joining method ("or" to include sources in any
                          catalogue, "and" to only include sources in all
                          catalogues)
        :type join_type: str
        :param mask: In preparation
        :type mask: list
        :param cartesian: if we should use cartesian coordinates
        :type cartesian: bool
        :param kwargs: extra args
        :type kwargs: dict
        :return: Matched catalogue
        :rtype: astropy.table.Table
        """
        catalogues = self.init_catalogues(catalogues)
        if self._col_names and TableColumn.CAT_NUM in self._col_names:
            self._col_names.remove(TableColumn.CAT_NUM)

        masked: Table = self.mask_catalogues(catalogues, mask)
        base: Table = self.build_meta(catalogues)

        if join_type == "and":
            p_error("join_type 'and' not fully implemented\n")

        # Bulk matching processes (column naming)
        for n, cat in enumerate(catalogues, 1):
            self._load.msg = "matching: %d" % n
            tmp: Table = self.inner_match(
                base, cat, join_type=join_type, cartesian=cartesian)

            tmp.rename_columns(
                tmp.colnames, [f"{name}_{n}" for name in tmp.colnames])
            # check if there is data to match against.
            if tmp is None or len(tmp) == 0:
                printf(f"No matches were found in catalogue {n}")
                continue
            base = fill_nan(hstack((base, tmp)))

        # Add in any masked bits
        if len(masked):
            masked.rename_columns(
                masked.colnames, [f"{name}_0" for name in masked.colnames])
            base = fill_nan(vstack((base, masked)))
        return base

    def inner_match(
            self, base: Table, cat: Table, join_type: str = "or",
            cartesian: bool = False) -> Table:
        """
        Base matching function between two catalogues

        :param base: The first catalogue for matching
        :type base: astropy.table.Table
        :param cat: The second catalogue for matching
        :type cat: astropy.table.Table
        :param join_type: Joining method ("or" to include all sources, "and"
                          for intersections)
        :type join_type: str
        :param cartesian: Whether to use Cartesian coordinates for matching
        :type cartesian: bool
        :return: Indices, 2D separation, and 3D separation
        :rtype: astropy.table.Table
        """
        if not len(base):
            return cat.copy()

        base = fill_nan(base.copy())
        assert self._col_names is not None
        col_names: list[str] = [
            n for n in self._col_names if n in cat.colnames]
        cat = fill_nan(cat[col_names].copy())

        sky_coords_1: SkyCoord
        sky_coords_2: SkyCoord
        if not cartesian:
            sky_coords_1 = self._sky_coords_not_cartesian(base)
            sky_coords_2 = self._sky_coords_not_cartesian(cat)
        else:
            sky_coords_1 = self._sky_coords_cartesian(base)
            sky_coords_2 = self._sky_coords_cartesian(cat)

        idx: np.ndarray
        d2d: units.Quantity
        d3d: units.Quantity
        idx, d2d, d3d = sky_coords_2.match_to_catalog_3d(sky_coords_1)

        tmp: Table = Table(
            np.full((len(base), len(col_names)), np.nan),
            names=col_names, dtype=cat[col_names].dtype)

        dist: np.ndarray | units.Quantity
        threshold: Any
        if cartesian:
            dist = d3d
            # If your threshold has an explicit unit wrapper,
            # extract its scalar magnitude
            threshold = self._threshold.value if hasattr(
                self._threshold, "value") else self._threshold
        else:
            dist = d2d
            threshold = self._threshold

        src: Any
        IDX: int
        sep: Any
        for src, IDX, sep in zip(cat.as_array(), idx, dist):
            self._load()
            if self._verbose:
                self._load.show()

            if (sep <= threshold) and (sep == min(dist[idx == IDX])):
                tmp[IDX] = src
            elif join_type == "or":
                tmp.add_row(src)

        return tmp

    def finish_matching(
            self,
            tab: Table | None,
            error_column: str = TableColumn.E_FLUX,
            num_thresh: int = -1,
            zp_mag: float = 0.0,
            col_names: list[str] | None = None) -> Table:
        """
        Averaging all the values. Combining source flags and building a NUM
        column

        :param tab: Table to work on
        :type tab: astropy.table.Table | None
        :param error_column: Column containing resultant photometric errors
                             (E_FLUX or STD_FLUX)
        :type error_column: str
        :param num_thresh: Minimum number of matches a source must have
                           (no cropping if <= 0)
        :type num_thresh: int
        :param zp_mag: Zero point (Magnitude) to be applied to the final
                       magnitude
        :type zp_mag: float
        :param col_names: List of column names to average; uses
                          self.col_names if None
        :type col_names: list, optional
        :return: An averaged version of the input table
        :rtype: astropy.table.Table
        """
        if tab is None:
            return Table(None)

        # have working table
        flags: np.ndarray = np.full(
            len(tab), SourceFlags.SRC_GOOD, dtype=np.uint16)
        average_table: Table = Table(None)

        if col_names is None:
            col_names = self._col_names if self._col_names else []

        for ii, name in enumerate(col_names):
            if all_cols := find_col_names(tab, name):
                ar: np.ndarray = tab2array(tab, col_names=all_cols)
                col: Column

                # only go forward if both tables have a common column name to
                # compare
                if ar.shape[1] > 1:
                    if name == TableColumn.FLUX:
                        col = Column(np.nanmedian(ar, axis=1), name=name)
                        mean: np.ndarray = np.nanmean(ar, axis=1)

                        if (self._col_names and
                                TableColumn.STD_FLUX not in self._col_names):
                            average_table.add_column(
                                Column(np.nanstd(ar, axis=1),
                                       name=TableColumn.STD_FLUX),
                                index=ii + 1)
                        # if median and mean are >5% different, flag as
                        # SRC_VAR.
                        # ABS. why are we using such an aggressive type check
                        # here?
                        flags[np.abs(mean - col) > (col / 5.0)] |= np.uint16(
                            SourceFlags.SRC_VAR)
                    elif name == TableColumn.E_FLUX:
                        col = Column(
                            np.sqrt(np.nansum(ar * ar, axis=1)), name=name)
                    elif name == TableColumn.STD_FLUX:
                        col = Column(np.nanmedian(ar, axis=1), name=name)
                    elif name == TableColumn.FLAG:
                        col = Column(flags, name=name)
                        f_col: np.ndarray
                        for f_col in ar.T:
                            flags |= f_col.astype(np.uint16)
                    elif name == TableColumn.NUM:
                        col = Column(np.nansum(ar, axis=1), name=name)
                    elif name == TableColumn.CAT_NUM:
                        col = Column(all_cols[0], name=name)
                    else:
                        col = Column(np.nanmedian(ar, axis=1), name=name)
                else:
                    col = tab[all_cols[0]]
                    col.name = name
                average_table.add_column(col, index=ii)

        average_table[TableColumn.FLAG] = Column(flags, name=TableColumn.FLAG)
        if TableColumn.FLUX in average_table.colnames:
            ecol: Column | None = (
                average_table[error_column]
                if error_column in average_table.colnames else None)
            mag: np.ndarray
            mag_err: np.ndarray
            mag, mag_err = flux_to_pogson_mag(
                average_table[TableColumn.FLUX], flux_err=ecol)
            mag += zp_mag

            if self._filter in average_table.colnames:
                average_table.remove_column(str(self._filter))
            if f"e{self._filter}" in average_table.colnames:
                average_table.remove_column(f"e{self._filter}")
            average_table.add_column(mag, name=str(self._filter))
            average_table.add_column(mag_err, name=f"e{self._filter}")

        if TableColumn.NUM not in average_table.colnames:
            narr: np.ndarray = np.nansum(np.invert(
                np.isnan(tab2array(
                    tab, find_col_names(tab, TableColumn.RA)))), axis=1)
            average_table.add_column(Column(narr, name=TableColumn.NUM))

            if num_thresh > 0:
                average_table.remove_rows(
                    average_table[TableColumn.NUM] < num_thresh)
        return average_table

    @property
    def col_names(self) -> list[str] | None:
        return self._col_names

    @property
    def filter(self) -> str | None:
        return self._filter

    @property
    def threshold(self) -> Any:
        return self._threshold

    @property
    def verbose(self) -> int | None:
        return self._verbose
