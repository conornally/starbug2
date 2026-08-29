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
import stat
import numpy as np
from typing import List, Optional, TextIO, Dict

from starbug2.constants import (
    FITS_EXTENSION, FILE_NAME, HeaderTags, ImageHeaderTags)
from astropy.io import fits
from starbug2.utilities.utils import printf, p_error, split_file_name

# A clear, Type Alias for the deep data nested structure Format:
# Dict[KeyType, ValueType], str mapping being FILTER, OBS, VISIT, DETECTOR
ExposureMapping = (
    Dict[
        Optional[int], Dict[
            Optional[int], Dict[
                Optional[int], Dict[
                    Optional[int], List[fits.HDUList]]]]])


def generate_runscript(
        f_names: List[str], args: str = "starbug2 ") -> None:
    """
    generate the run script

    :param f_names: file names
    :type f_names: list of str.
    :param args: arguments
    :type args: str
    :return: None
    """
    runfile: str = "./run.sh"
    fits_files: List[fits.HDUList] = []

    fp: TextIO = open(runfile, "w")
    fp.write("#!/bin/bash\n")
    fp.write("CMDS=\"-vf\"\n")
    for f_name in f_names:
        if os.path.exists(f_name):
            d_name: str
            name: str
            ext: str
            d_name, name, ext = split_file_name(f_name)
            if ext == FITS_EXTENSION:
                fits_file: fits.HDUList = fits.open(f_name)
                fits_file[0].header[FILE_NAME] = f_name
                fits_files.append(fits_file)
            else:
                p_error(
                    "file %s must be type '.fits' not '%s'\n" % (name, ext))
        else:
            p_error("file \x1b[1;31m%s\x1b[0m not found\n" % f_name)

    sorted_exposures: ExposureMapping = sort_exposures(fits_files)

    # print exps.
    for _, obs in sorted_exposures.items():
        for _, visits in obs.items():
            for _, destinations in visits.items():
                for _, exps in destinations.items():
                    s = f"{args}${{CMDS}} -n{len(exps)} "
                    for exp in exps:
                        s += "%s " % exp[0].header[FILE_NAME]
                    fp.write("%s\n" % s)
    fp.close()
    os.chmod(
        runfile, stat.S_IXUSR | stat.S_IWUSR | stat.S_IRUSR | stat.S_IRGRP |
        stat.S_IROTH)
    printf("->%s\n" % runfile)


def sort_exposures(catalogues: List[fits.HDUList]) -> ExposureMapping:
    """
     Given a list of catalogue files, this will return the fitsHDULists as a
     series of nested dictionaries sorted by:
    >   BAND
    >   OBSERVATION ID
    >   VISIT ID
    >   DETECTOR                -- These two have been switched
    >   DITHER (EXPOSURE)       -- These two have been switched

    :param catalogues: the catalogues to sort exposures of.
    :type catalogues: list of fits.HDUList
    :return: a dictionary of sorted catalogues
    :rtype: ExposureMapping
    """
    out: ExposureMapping = {}
    for cat in catalogues:
        info = exp_info(cat)

        if info[HeaderTags.FILTER] not in out.keys():
            out[info[HeaderTags.FILTER]] = {}

        if info[HeaderTags.OBS] not in out[info[HeaderTags.FILTER]].keys():
            out[info[HeaderTags.FILTER]][info[HeaderTags.OBS]] = {}

        if (info[HeaderTags.VISIT] not in
                out[info[HeaderTags.FILTER]][info[HeaderTags.OBS]].keys()):
            out[info[HeaderTags.FILTER]][
                info[HeaderTags.OBS]][info[HeaderTags.VISIT]] = {}

        if (info[ImageHeaderTags.DETECTOR] not in
            out[info[HeaderTags.FILTER]][info[HeaderTags.OBS]][
                info[HeaderTags.VISIT]].keys()):
            out[info[HeaderTags.FILTER]][
                info[HeaderTags.OBS]][info[HeaderTags.VISIT]][
                info[ImageHeaderTags.DETECTOR]] = []
        out[info[HeaderTags.FILTER]][
            info[HeaderTags.OBS]][info[HeaderTags.VISIT]][
            info[ImageHeaderTags.DETECTOR]].append(cat)
    return out


def parse_mask(string, table) -> np.ndarray | None:
    """
    Parse a commandline mask string to be passed into a matching routine
    Example: --mask=F444W!=nan

    :param string: Raw mask sting to be parsed
    :type string: str
    :param table: Table to work on
    :type table: astropy.table.Table
    :return: Boolean mask array to index into a table or array
    :rtype: np.ndarray
    """
    mask: Optional[np.ndarray] = None

    col_name: str
    for col_name in table.colnames:
        string: str = string.replace(col_name, "table[\"%s\"]" % col_name)

    try:
        mask = eval(string)
        if not isinstance(mask, np.ndarray):
            raise Exception
    except NameError as e:
        p_error("Unable to create mask: %s\n" % repr(e))
    except Exception as e:
        p_error(repr(e))

    return mask


def exp_info(hdu_list) -> Dict[str, int | None]:
    """
    Get the exposure information about a hdu list
    :param hdu_list: HDUList or ImageHDU or BinTableHDU
    :return: dictionary of relevant information
    (HeaderTags.FILTER, obs, visit exposure, detector)
    :rtype dict(str, Optional[int])
    """
    info: Dict[str, int | None] = {
        HeaderTags.FILTER: None,
        HeaderTags.OBS: 0,
        HeaderTags.VISIT: 0,
        HeaderTags.EXPOSURE: 0,
        ImageHeaderTags.DETECTOR: None
    }

    if type(hdu_list) in (fits.ImageHDU, fits.BinTableHDU):
        hdu_list: fits.HDUList = fits.HDUList(hdu_list)

    for hdu in hdu_list:
        for key in info:
            if key in hdu.header:
                info[key] = hdu.header[key]
    return info
