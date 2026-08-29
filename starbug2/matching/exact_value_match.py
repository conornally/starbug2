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
from typing import override, Any

import numpy as np
from astropy.table import Table, hstack, vstack

from starbug2.constants import TableColumn
from starbug2.matching.generic_match import GenericMatch
from starbug2.utilities.utils import p_error, fill_nan


class ExactValueMatch(GenericMatch):
    """
    Match catalogues based on a single value in one of their columns
    The normal use case is matching sources based on their *Catalogue_Number*
    value.
    """

    def __init__(
            self, value: str = TableColumn.CAT_NUM, **kwargs: Any) -> None:
        """
        Setup method.

        :param value: Column name to take exact values from
        :type value: str
        :param kwargs:
        """
        self.value: str = value
        super().__init__(**kwargs, method="Exact Value Matching")

        # noinspection SpellCheckingInspection
        if "colnames" in kwargs:
            # noinspection SpellCheckingInspection
            p_error("Colnames not implemented in %s\n" % self.method)

    def __str__(self) -> str:
        # noinspection SpellCheckingInspection
        s: list[str] = [
            f"{self.method}:",
            f'Value: "{self.value}"',
            f"Colnames: {self._col_names}",
        ]
        return "\n".join(s)

    @override
    def inner_match(
            self, base: Table, cat: Table, join_type: str = "or",
            cartesian: bool = False) -> Table:
        """
         The low level matching function.

        :param base: Table onto which to match *cat*
        :type base: astropy.Table
        :param cat: the Table to match to *base*
        :type cat:  astropy.Table
        :param join_type: Joining method ("or" to include all sources, "and"
                          for intersections)
        :type join_type: str
        :param cartesian: Whether to use Cartesian coordinates for matching
        :type cartesian: bool.
        :return: A new catalogue, it is a reordered version of *cat*, in the
            correct sorting to be h-stacked with *base*.
        :rtype:  astropy.Table.
        """
        tmp: Table = Table(
            np.full((len(base), len(cat.colnames)), np.nan),
            names=cat.colnames, dtype=cat.dtype, masked=True
        )

        for col in tmp.columns.values():
            col.mask = True

        if not len(base):
            return vstack([tmp, cat])

        for src in cat:
            if self._verbose:
                self._load()
                self._load.show()

            ii: np.ndarray = np.where(base[self.value] == src[self.value])[0]
            if len(ii):
                tmp[ii] = src
            else:
                tmp.add_row(src)
        return tmp

    @override
    def match(self, catalogues: list[Table], **kwargs: Any) -> Table | None:
        """
        Core matching function

        :param catalogues: List of astropy Tables to match together
        :type catalogues: list of astropy.Table
        :param kwargs:
        :return: Full matched catalogue
        :rtype: astropy.Table
        """
        catalogues = self.init_catalogues(catalogues)
        base: Table = self.build_meta(catalogues)

        if self._col_names is None or self.value not in self._col_names:
            p_error("Exact value '%s' not in column names.\n" % self.value)
            return None

        for n, cat in enumerate(catalogues, 1):
            self._load.msg = "matching: %d" % n
            tmp: Table = self.inner_match(base, cat)
            tmp.rename_columns(
                tmp.colnames, [f"{name}_{n}" for name in tmp.colnames]
            )
            base = hstack([base, tmp])

            if n > 1:
                ii: np.ndarray = (
                    base[self.value].mask & ~base[f"{self.value}_{n}"].mask
                )
                base[self.value][ii] = base[f"{self.value}_{n}"][ii]
                base.remove_column(f"{self.value}_{n}")
            else:
                base.rename_column(f"{self.value}_1", self.value)

        return fill_nan(base)
