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
from astropy.table import Table
from starbug2.constants import TableColumn
from starbug2.matching.generic_match import GenericMatch
from starbug2.utilities.utils import h_cascade, fill_nan


class CascadeMatch(GenericMatch):
    """
    A simple advancement on "Generic Matching" where the number
    of columns are not preserved. At the end of each sub match, the
    table is left justified, to reduce the total number of columns needed.
    """

    def __init__(self, **kwargs: Any) -> None:
        super().__init__(**kwargs, method="Cascade Matching")

    @override
    def match(self, catalogues: list[Table], **kwargs: Any) -> Table:
        """
        Match a list of catalogues with RA and DEC columns.

        :param catalogues: The input catalogues to work on
        :type catalogues: list (astropy.Tables)
        :param kwargs: Keyword arguments passed to GenericMatch.match
        :return: A left aligned catalogue of all the matched values
        :rtype: astropy.table.Table
        """
        catalogues = self.init_catalogues(catalogues)
        if self._col_names and TableColumn.CAT_NUM in self._col_names:
            self._col_names.remove(TableColumn.CAT_NUM)

        base: Table = self.build_meta(catalogues)

        for n, cat in enumerate(catalogues, 1):
            self._load.msg = "matching: %d" % n
            tmp: Table = self.inner_match(base, cat, join_type="or")
            tmp.rename_columns(
                tmp.colnames, [f"{name}_{n}" for name in tmp.colnames]
            )
            base = h_cascade([base, tmp], col_names=self._col_names)

        base = fill_nan(base)
        return base
