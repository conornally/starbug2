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
import sys
import time
from importlib import metadata
from importlib.metadata import PackageNotFoundError
from typing import Any, Dict, List, Tuple

from astropy.io import fits
from astropy.table import Column, MaskedColumn, Table, hstack, vstack
from astropy.wcs import WCS
import numpy as np
import requests

from starbug2.constants import (
    DEFAULT_COLOUR,
    TMP_OUT,
    TMP_FITS,
    TableColumn,
    HeaderTags,
    FITS_EXTENSION,
    N_MIS_MATCHES,
    ExitStates,
    REST_SUCCESS_CODE,
    Units,
    ImageHeaderTags, STARBUG_DATA_DIR, DEFAULT_BUNIT,
)
from starbug2.utilities.filters import STAR_BUG_FILTERS


# different print methods (why are we not using loggers?)
def printf(s: str) -> int:
    if s.endswith("\n"):
        return sys.stdout.write(s)
    else:
        return puts(s)


def p_error(s: str) -> int:
    return sys.stderr.write(s)


def puts(s: str) -> int:
    return printf("%s\n" % s)


def s_bold(s: str) -> str:
    return "\x1b[1m%s\x1b[0m" % s


def warn(s: str) -> int:
    return p_error("%s%s" % (s_bold("Warning: "), s))


def append_chars(s: str, n: int, c: str) -> str:
    """Append n characters to s.

    :param s: the base string
    :type s: str
    :param n: the number of times to add the character
    :type n: int
    :param c: the characters to add.
    :type c: str
    :return: the adjusted string.
    :rtype: str
    """
    for _ in range(n):
        s += c
    return s


def repeat_print(n: int, c: str) -> None:
    """Prints out a repeated string.

    NOTE: this seems unused.

    :param n: the number of times to repeat
    :type n: int
    :param c: the string to repeat.
    :type c: str
    :return: None
    """
    printf(append_chars("", n, c))


def split_file_name(file_path: str | None) -> Tuple[str, str, str]:
    """Breaks apart a path into folder, filename and extension.

    :param file_path: the path to split
    :return: (folder, file name, extension)
    :rtype: tuple of str, str, str
    """
    if file_path is None:
        raise Exception("failed as path is None")

    folder, file = os.path.split(file_path)
    file_name, ext = os.path.splitext(file)
    if not folder:
        folder = "."
    return folder, file_name, ext


class Loading(object):
    # how long the bar is
    bar = 40

    # current length
    n = 0

    # no idea
    length = 1

    # loading bar message
    msg = ""

    def __init__(self, length: int, msg: str = "", res: int = 1) -> None:
        self.set_len(length)
        self.msg = msg
        self.start_time = time.time()
        self.res = res

    def set_len(self, length: int) -> None:
        self.length = abs(length)

    def __call__(self) -> bool:
        self.n += 1
        return self.n <= self.length

    def show(self) -> None:
        dec: int = int(self.n / self.length)
        # only show once per self.res loads
        if (dec == 1) or (not self.n % self.res):
            out: str = "%s|" % self.msg
            for i in range(self.bar + 0):
                out += "=" if (i < (self.bar * dec)) else " "
            out += "|%.0f%%" % (100 * dec)

            if self.n:
                etc: float = (
                    (time.time() - self.start_time)
                    * (self.length - self.n)
                    / self.n
                )
                n_hrs: float = etc // 3600
                n_minutes: float = (etc - (n_hrs * 3600)) // 60
                n_secs: float = etc - (n_hrs * 3600) - (n_minutes * 60)
                stime: str = ""
                if n_hrs:
                    stime += "%dh" % int(n_hrs)
                if n_minutes:
                    stime += "%dm" % int(n_minutes)

                stime += "%ds" % int(n_secs)
                out += " ETC:%s" % stime

            printf("\x1b[2K%s\r" % out)
            sys.stdout.flush()
        if dec == 1:
            printf("\n")


def combine_tables(base: Table | None, tab: Table | None) -> Table | None:
    """Is this the same as vstack?"""
    if not base:
        return tab
    else:
        return vstack([base, tab])


def export_region(
    tab: Table,
    colour: str = DEFAULT_COLOUR,
    scale_radius: int = 1,
    region_radius: int = 3,
    x_col: str = TableColumn.RA,
    y_col: str = TableColumn.DEC,
    wcs: int = 1,
    f_name: str = TMP_OUT,
) -> None:
    """A handy function to convert the detections in a DS9 region file.

    :param tab: Source list table with some kind of positional columns
    :type tab: Table
    :param colour: Region Colour
    :type colour: str
    :param scale_radius: Scale region radius with flux ? true/false
    :type scale_radius: int
    :param region_radius: Otherwise, use this region radius in pixels
    :type region_radius: int
    :param x_col: X column name to use
    :type x_col: str
    :param y_col: Y column name to use
    :type y_col: str.
    :param wcs: Boolean which is true if the x or y_cols use the WCS system.
    :type wcs: int.
    :param f_name: Filename to output to
    :type f_name: str
    :return:
    """
    if x_col not in tab.colnames:
        x_cols = list(filter(lambda s: "x" == s[0], tab.colnames))
        if x_cols:
            x_col = x_cols[0]
            printf("Using '%s' as x position column\n" % s_bold(x_col))
            wcs = 0

    if y_col not in tab.colnames:
        y_cols = list(filter(lambda s: "y" == s[0], tab.colnames))
        if y_cols:
            y_col = y_cols[0]
            printf("Using '%s' as y position column\n" % s_bold(y_col))
            wcs = 0

    r: np.ndarray
    if TableColumn.FLUX in tab.colnames and scale_radius:
        r = -40.0 / np.log10(tab[TableColumn.FLUX])
        r[r < region_radius] = region_radius
        r[np.isnan(r)] = region_radius
    else:
        r = np.ones(len(tab)) * region_radius

    prefix: str = "fk5;" if wcs else ""

    with open(f_name, "w") as fp:
        fp.write("global color=%s width=2\n" % colour)
        if tab:
            for src, ri in zip([tab], r[r > 0]):
                fp.write(
                    "%scircle %f %f %fi\n" % (
                        prefix, src[x_col], src[y_col], ri)
                )
        else:
            p_error("unable to open %f\n" % f_name)


def translate_param_float(
    opt: str,
    opt_arg: str,
    set_opt: Dict[str, float],
    options: int,
    kill_option: int,
) -> Tuple[int, Dict[str, float]]:
    """Converts an opt param into a float.

    :param opt: the opt string
    :type opt: str
    :param opt_arg: the opt_arg string to convert
    :type opt_arg: str
    :param set_opt: the set opt dictionary
    :type set_opt: dict of strings
    :param options: the options integer
    :type options: int
    :param kill_option: the kill bit mask
    :type kill_option: int
    :return:  tuple of options int and set_opt dict
    :rtype int, dict
    """
    if opt in ("-s", "--set"):
        if "=" in opt_arg:
            key: str
            val: float
            key, val = opt_arg.split("=")
            try:
                val = float(val)
            except ValueError:
                pass
            set_opt[key] = val
        else:
            p_error("unable to set parameter, use syntax -s KEY=VALUE\n")
            options |= kill_option
    return options, set_opt


def parse_unit(raw: str) -> Tuple[float | None, int | None]:
    """Take a value with the ability to be cast into several units and
    parse it.

    i.e. 123p -> 123 'pixels'

    Recognised units are:
    p : pixels
    s : arcsec
    m : arcmin
    d : degree

    :param raw: Raw input string to operate on
    :type raw: str
    :return: Numerical value of unit
    :rtype float
    """
    recognised: Dict[str, int] = {
        "p": Units.PIX,
        "s": Units.ARCSEC,
        "m": Units.ARCMIN,
        "d": Units.DEG,
    }
    value: float | None = None
    unit: int | None = None
    if raw:
        try:
            value = float(raw)
            unit = None
        except ValueError:
            try:
                value = float(raw[:-1])
                unit = recognised.get(raw[-1])
            except ValueError:
                p_error("unable to parse '%s'\n" % raw)
    return value, unit


def tab2array(tab: Table, col_names: List[str] | None = None) -> np.ndarray:
    """Returns the contents of the table as a normal 2D numpy array.

    NB: this is different from Table.asarray(), which returns an array of
    numpy.voids

    if col_names not None, return the subset of the table corresponding to
    this list

    :param tab: Table to operate on
    :type tab: astropy.Table
    :param col_names: Column names in table to include in the array
    :type col_names: list of str
    :return: Array from the table
    :rtype: numpy.ndarray
    """
    if not col_names:
        col_names = tab.colnames
    else:
        col_names = remove_duplicates(col_names)
    return np.array(tab[col_names].as_array().tolist())


def collapse_header(header) -> fits.Header:
    """Convert a dictionary to a Header.

    Parameters in PARAMFILES have keys longer than 8 chars
    which can cause issues in the fits format. This function turns
    those to comment cards.

    :param header: Header or dictionary to convert to collapse header
    :type header: dict, fits.Header
    :return: Collapsed Header
    :rtype fits.Header
    """
    out: fits.Header = fits.Header()
    key: str
    value: float
    for key, value in header.items():
        if len(key) > 8:
            out["comment"] = ":".join([key, str(value)])
        else:
            out[key] = value
    return out


def export_table(
    table: Table,
    f_name: str | None = None,
    header: fits.Header | None = None,
) -> None:
    """Export table with correct dtypes.

    :param table: Table to export.
    :type table: astropy.Table
    :param f_name: Filename to export to.
    :type f_name: str
    :param header: Optional header file to include in fits table
    :type header: dict, fits.Header
    :return: None
    """
    dtypes: List[Any] = []
    if TableColumn.CAT_NUM not in table.colnames:
        table = reindex(table)
    for name in table.colnames:
        if name == TableColumn.CAT_NUM:
            dtypes.append(str)
        elif name == "flag":
            dtypes.append(np.uint16)
        else:
            dtypes.append(table[name].dtype)
    table = fill_nan(Table(table, dtype=dtypes))

    if not f_name:
        f_name = TMP_FITS

    fits_header = header if header is not None else fits.Header()

    fits.BinTableHDU(data=table, header=fits_header).writeto(
        f_name, overwrite=True, output_verify="fix"
    )


def import_table(f_name: str, verbose: bool | int = 0) -> Table | None:
    """Slight tweak to `astropy.table.Table.read`.

    This makes sure that the proper column dtypes are maintained

    :param f_name: Path to binary fits table file
    :type f_name: str
    :param verbose: Display verbose information
    :type verbose: boolean or int
    :return: Loading table
    :rtype: atrophy.Table | None
    """
    printf(f"trying to load file {f_name}\n")
    printf(f" file exists in correct location {os.path.exists(f_name)}\n")
    new_f_name = os.path.join(
        str(os.environ.get("STARBUG_DATDIR")), os.path.basename(f_name))
    printf(f"new location path is {new_f_name}\n")
    printf(f"the file exists in datadir instead {
        os.path.exists(new_f_name)}\n")

    tab: Table | None = None
    if os.path.exists(f_name):
        if os.path.splitext(f_name)[1] == FITS_EXTENSION:
            tab = fill_nan(Table.read(f_name, format="fits"))
            if tab is None:
                printf(f"table at {f_name} failed to read")
                return None
            if not tab.meta.get(HeaderTags.FILTER):
                if filter_string := find_filter(tab):
                    tab.meta[HeaderTags.FILTER] = filter_string
            if verbose:
                printf(
                    "-> loaded %s (%s:%d)\n"
                    % (f_name, tab.meta.get(HeaderTags.FILTER), len(tab))
                )
        else:
            p_error("Table must fits format\n")
    else:
        p_error("Unable to locate \"%s\"\n" % f_name)
    return tab


def fill_nan(table: Table) -> Table:
    """Fill empty values in table with nans.

    This is useful for tables that have columns that don't support nans (e.g.
    starbug flag). These will be set to zero instead

    :param table: table to operate on
    :type table: atrophy.table
    :return: Input table with masked vales filled in as nan
    :rtype: atrophy.table
    """
    name: str
    fill_val: int | float
    for _, name in enumerate(table.colnames):
        match table[name].dtype.kind:
            case "f":
                fill_val = np.nan
            case "i" | "u":
                fill_val = 0
            case _:
                fill_val = np.nan
        if isinstance(table[name], MaskedColumn):
            table[name] = table[name].filled(fill_val)
    return table


def find_col_names(tab: Table, basename: str) -> List[str]:
    """Find substring (basename) within the table colnames.

    Searches for substring at the beginning of the word I.E search for "flux"
    in ("flux_out","flux_err","d_flux") returns as ("flux_out","flux_err")

    :param tab: Table to operate on
    :type tab: atrophy.table
    :param basename: String basename to search
    :type basename: str
    :return: List of all matching column names
    :rtype: list of str
    """
    return [
        col_name
        for col_name in tab.colnames
        if col_name[: len(basename)] == basename
    ]


def combine_file_names(
    f_names: List[str | None], n_mismatch: int = N_MIS_MATCHES
) -> str | None:
    """When matching catalogues, combines the file names.

    Combines file names into an appropriate combination of all the inputs.

    :param f_names: list of file names
    :type f_names: list of str | None
    :param n_mismatch: The number of mismatched characters it will allow
    :type n_mismatch: int
    :return: Combined filenames
    :rtype: str
    """
    trys: int = 0
    f_name: str = ""
    d_name: str
    ext: str

    if f_names is None:
        return None

    d_name, _, ext = split_file_name(f_names[0])
    f_names_split: List[str] = [split_file_name(name)[1] for name in f_names]

    for i in range(len(f_names_split[0])):
        chars: List[str] = [
            name[i] for name in f_names_split if len(name) > i
        ]
        if len(set(chars)) == 1:
            f_name += chars[0]
        else:
            f_name += "(%s)" % "".join(sorted(set(chars)))
            trys += 1
        if trys > n_mismatch:
            return None
    while ")(" in f_name:
        f_name = f_name.replace(")(", "")
    return "%s/%s%s" % (d_name, f_name, ext)


def h_cascade(
    tables: List[Table], col_names: List[str] | None = None
) -> Table:
    """Similar use as hstack.

    Except rather than adding a full new column, the inserted value is placed
    into the leftmost empty column

    :param tables: Table to h_cascade.
    :type tables: list of atrophy.Table
    :param col_names: List of column names to include in the stacking.
        If col_names=None, use all possible columns
    :type col_names: list of str
    :return: Single combined table
    :rtype: atrophy.Table
    """
    tab: Table = fill_nan(hstack(tables))

    if not col_names:
        col_names = tables[0].colnames
    for name in col_names:
        cols: List[str] = find_col_names(tab, name)
        if not cols:
            continue
        move: int = 1
        while move:
            move = 0
            for n in range(len(cols) - 1, 0, -1):
                curr_mask: np.ndarray = np.invert(np.isnan(tab[cols[n]]))
                left_mask: np.ndarray = np.isnan(tab[cols[n - 1]])
                mask: np.ndarray = np.logical_and(curr_mask, left_mask)

                tab[cols[n]] = MaskedColumn(tab[cols[n]])
                tab[cols[n - 1]][mask] = tab[cols[n]][mask]
                tab[cols[n]][mask] = tab[cols[n]].info.mask_val
                if sum(mask):
                    move = 1
        cols = find_col_names(tab, name)
        if cols:
            tab.rename_columns(
                cols, ["%s_%d" % (name, i + 1) for i in range(len(cols))]
            )

    for name in tab.colnames:
        col: Table = tab[name]
        n_bad: int = getattr(col.info, "n_bad", 0)
        if n_bad == len(col):
            tab.remove_column(name)
    return tab


def ext_names(hdu_list: fits.HDUList | None) -> List[str]:
    """Return list of HDU extension names.

    :param hdu_list: fits hdu_list to operate on
    :type hdu_list: fits.HDUList
    :return: List of extension names
    :rtype: list of str
    """
    if hdu_list is None:
        return []

    return [ext.name for ext in hdu_list]


# noinspection SpellCheckingInspection
def flux_to_pogson_mag(
    raw_flux: np.ndarray | float,
    flux_err: Column | None | np.ndarray | int | float = None
) -> Tuple[np.ndarray, np.ndarray]:
    """Convert flux to magnitude in an arbitrary system.

    Uses the Pogsons relation.

    :param raw_flux: List of source flux values
    :type raw_flux: list of floats or float or None or ndarray
    :param flux_err: List of known flux uncertainties
    :type flux_err: list of floats or float or None or ndarray
    :return: tuple of (Source magnitudes, Magnitude errors ).
    :rtype: tuple (ndarray, ndarray).
    """
    flux: np.ndarray = np.atleast_1d(raw_flux).astype(float)
    if flux_err is None:
        flux_err_arr: np.ndarray = np.zeros_like(flux)
    else:
        flux_err_arr = np.atleast_1d(np.asarray(flux_err)).astype(float)
    mag: np.ndarray = np.full(len(flux), np.nan)
    mag_err: np.ndarray = np.full(len(flux), np.nan)

    with np.errstate(invalid="ignore"):
        mask_flux: np.ndarray = (flux > 0) & np.isfinite(flux)
        mask_f_err: np.ndarray = flux_err_arr >= 0

        mask: np.ndarray = mask_flux & mask_f_err

        mag[mask_flux] = -2.5 * np.log10(flux[mask_flux])
        mag_err[mask] = 2.5 * np.log10(1.0 + (flux_err_arr[mask] / flux[mask]))
    mag[flux == np.inf] = -np.inf

    return mag, mag_err


def wget(address: str, f_name: str | None = None) -> ExitStates:
    """A really simple "implementation" of wget.

    :param address: URL to download
    :type address: str
    :param f_name: Filename to save output to
    :type f_name: str
    :return: 0 on success, 1 on failure
    :rtype ExitStates
    """
    r: requests.Response = requests.get(address)
    if r.status_code == REST_SUCCESS_CODE:
        f_name = f_name if f_name else os.path.basename(address)
        with open(f_name, "wb") as fp:
            for chunk in r.iter_content(chunk_size=128):
                fp.write(chunk)
        return ExitStates.EXIT_SUCCESS
    else:
        p_error("Unable to download \"%s\"\n" % address)
        return ExitStates.EXIT_FAIL


def reindex(table: Table) -> Table:
    """Add indexes into a table.

    :param table: the table to reindex
    :type table: atrophy.Table
    :return: the reindex-ed table
    :rtype: atrophy.Table
    """
    if TableColumn.CAT_NUM in table.colnames:
        table.remove_column(TableColumn.CAT_NUM)
    column = Column(
        ["CN%d" % i for i in range(len(table))], name=TableColumn.CAT_NUM
    )
    table.add_column(column, index=0)
    return table


def colour_index(table: Table, keys: List[str]) -> Table:
    """Allow table indexing with A-B.

    :param table: table to colour index.
    :type table: atrophy.Table
    :param keys: column names to index.
    :return: A table which has only columns defined in keys.
    :rtype: atrophy.Table
    """
    out: Table = Table()
    key: str
    for key in keys:
        if key in table.colnames:
            out.add_column(table[key])
        elif "-" in key:
            a: str
            b: str
            a, b = key.split("-")
            out.add_column(table[a] - table[b], name=key)
    return out


def get_mj_ysr2jy_scale_factor(
    ext: fits.PrimaryHDU | fits.ImageHDU | fits.BinTableHDU,
) -> float:
    """Find the unit scale factor to convert an image from MJy/sr to Jy.

    Header file must contain the keyword "PIXAR_SR"

    :param ext: Fits extension with header file
    :type ext: PrimaryHDU, ImageHDU, BinTableHDU
    :return: Value of scaling factor from the header
    :rtype float
    """
    scale_factor: float = 1.0
    if ext.header.get(ImageHeaderTags.BUN_IT) == DEFAULT_BUNIT:
        if ImageHeaderTags.PIXAR_SR in ext.header:
            scale_factor = 1e6 * float(ext.header[ImageHeaderTags.PIXAR_SR])
    return scale_factor


def find_filter(table: Table) -> str:
    """Attempt to identify filter for a table from metadata or columns.

    :param table: Table to work on.
    :type table: astropy.table.Table
    :return: Identified filter value, otherwise None.
    :rtype: str
    """
    filter_string: str
    if filter_string := table.meta.get(HeaderTags.FILTER):
        return filter_string

    matching_filters: set[str] = (
        set(table.colnames) & set(STAR_BUG_FILTERS.keys())
    )
    if matching_filters:
        return matching_filters.pop()

    return ""


def get_version() -> str:
    """Try to determine the installed starbug version on the system.

    :return: the StarBugII version string
    :rtype str
    """
    version: str
    try:
        version = metadata.version("starbug2")
    except (AttributeError, TypeError, PackageNotFoundError):
        version = "UNKNOWN"
    return version


def remove_duplicates[T](seq: List[T]) -> List[T]:
    """Take a sequence and rm its duplicates while preserving order.

    :param seq: Input list to work on
    :type seq: list of <T>
    :return: A copy of the list with the duplicate elements removed
    :rtype list of <T>
    """
    seen: set[T] = set()
    seen.update(seq)
    to_return: list[T] = []
    value: T
    for value in seq:
        if value in seen:
            to_return.append(value)
            seen.remove(value)
    return to_return


def crop_hdu(
    hdu: fits.HDUList,
    x_limit: Tuple[int, int] | None = None,
    y_limit: Tuple[int, int] | None = None,
) -> fits.HDUList | None:
    """Crop an image with multiple extensions. Retaining the extensions.

    :param hdu: A multi frame fits HDUList
    :type hdu: fits.HDUList
    :param x_limit: Pixel X bounds to crop image between
    :type x_limit: list of ?????
    :param y_limit: Pixel Y bounds to crop image
    :type y_limit: list of ????
    :return: The full HDUList that has been spatially cropped
    :rtype fits.HDUList
    """
    if x_limit is None or y_limit is None:
        return None

    ext: fits.PrimaryHDU | fits.ImageHDU | fits.BinTableHDU | fits.GroupsHDU
    for ext in hdu:
        if not isinstance(ext, (fits.PrimaryHDU, fits.ImageHDU)):
            continue
        if not ext.header[HeaderTags.NAXIS]:
            continue

        ctype: str = ext.header.get(HeaderTags.C_TYPE)
        ext.header[HeaderTags.C_TYPE] = "%s-SIP" % ctype

        w: WCS = WCS(ext.header, relax=False)
        ext.data = ext.data[x_limit[0]: x_limit[1], y_limit[0]: y_limit[1]]
        ext.header.update(
            w[x_limit[0]: x_limit[1], y_limit[0]: y_limit[1]].to_header()
        )
    return hdu


def usage(docstring: str | None, verbose: bool | int = 0) -> int:
    """Outputs the usage layout string.

    :param docstring: the doc string to output
    :param verbose: if to do so in verbose mode
    :return: 1 when complete.
    """
    if docstring is None:
        return 1

    if verbose:
        p_error(docstring)
    else:
        p_error("%s\n" % docstring.split("\n")[1])
    return 1


def parse_cmd(args: List[str]) -> Tuple[str, List[str]]:
    """Parses an args command.

    :param args: the args array.
    :return: tuple of the command and the rest of the args array.
    :rtype: (str, array[str])
    """
    cmd = os.path.basename(args[0])
    return cmd, args[1:]


def get_data_path() -> str:
    """
    Returns the data path.

    :return: The data path
    :rtype: str
    """
    env_path: str | None = os.getenv(STARBUG_DATA_DIR)
    return (env_path if env_path else
            "%s/.local/share/starbug" % (os.getenv("HOME")))


if __name__ == "__main__":
    print(parse_unit(""))
    print(parse_unit("10p"))
    print(parse_unit("10 p"))
    print(parse_unit("10 D"))
    print(parse_unit("10"))
    print(parse_unit("p10"))
