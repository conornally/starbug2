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

from typing import List, Final
from enum import Enum
from pathlib import Path
import os

# the filter id which we've had to adjust the bin size to allow it to
# initialise without errors.
PROBLEMATIC_FILTER_ID = "F150W2"
PROBLEMATIC_FILTER_WARNING = (
    "Caution needed with F150W2 photometric accuracy. More info"
    "can be found in ("
    "https://github.com/alan-stokes/starbug2/issues/2)")

# noinspection SpellCheckingInspection
STARBUG_DATA_DIR: Final[str] = "STARBUG_DATDIR"
# noinspection SpellCheckingInspection
WEBBPSF_PATH_ENV_VAR: Final[str] = "WEBBPSF_PATH"
# noinspection SpellCheckingInspection
STAR_BUG_PARAMS: Final[str] = "STARBUGII PARAMETERS"
STAR_BUG_TEST_DAT_ENV: Final[str] = "STARBUG_TEST_DIR"

# default values
DEFAULT_FULL_WIDTH_HALF_MAX = 2.0
DEFAULT_PSF_FILE_NAME = "psf.fits"
DEFAULT_COLOUR: Final[str] = "green"
DEFAULT_MIN_MAG: Final[int] = 28
DEFAULT_MAX_MAG: Final[int] = 18
DEFAULT_BUNIT = "MJy/sr"

# how many characters we will allow by default.
N_MIS_MATCHES: Final[int] = 10

# rest success
REST_SUCCESS_CODE: Final[int] = 200

# url to docs
URL_DOCS: Final[str] = (
    "https://raw.githubusercontent.com/conornally/starbug2/"
    "refs/heads/main/docs/source/_static/images/starbug.png")

READ_THE_DOCS_URL: Final[str] = "https://starbug2.readthedocs.io/en/latest/"

# fit urls
# noinspection SpellCheckingInspection
JWST_MIRI_APCORR_0010_FITS_URL: Final[str] = (
    "https://jwst-crds.stsci.edu/unchecked_get/references/jwst/"
    "jwst_miri_apcorr_0010.fits"
)
# noinspection SpellCheckingInspection
JWST_NIRCAM_APCORR_0004_FITS_URL: Final[str] = (
    "https://jwst-crds.stsci.edu/unchecked_get/references/jwst/"
    "jwst_nircam_apcorr_0004.fits"
)
# noinspection SpellCheckingInspection
# ab vega offset urls
JWST_MIRI_ABVEGA_OFFSET_URL: Final[str] = (
    "https://jwst-crds.stsci.edu/unchecked_get/references/jwst/"
    "jwst_miri_abvegaoffset_0001.asdf"
)
# noinspection SpellCheckingInspection
JWST_NIRCAM_ABVEGA_OFFSET_URL: Final[str] = (
    "https://jwst-crds.stsci.edu/unchecked_get/references/jwst/"
    "jwst_nircam_abvegaoffset_0002.asdf"
)

# paths to temp files.
TMP_OUT: Final[str] = "/tmp/out.reg"
TMP_FITS: Final[str] = "/tmp/starbug.fits"

# the fits file extension
FITS_EXTENSION: Final[str] = ".fits"
FILE_NAME: Final[str] = "FILENAME"

# HDU extension names
DQ: Final[str] = "DQ"
AREA: Final[str] = "AREA"
WHT: Final[str] = "WHT"
ERR: Final[str] = "ERR"

# file types
AP_FILE: Final[str] = "AP_FILE"
BGD_FILE: Final[str] = "BGD_FILE"
PSF_FILE: Final[str] = "PSF_FILE"

# the number of columns in the test table.
N_COLUMNS: Final[int] = 8

# flags inside the artificial stars results table.
NOT_FOUND: Final[int] = 0
DETECT: Final[int] = 1


# SOURCE FLAGS
class SourceFlags(int, Enum):
    SRC_GOOD = 0
    SRC_BAD = 0x01
    SRC_JMP = 0x02
    # source frame mean >5% different from median
    SRC_VAR = 0x04
    # psf fit with fixed centroid
    SRC_FIX = 0x08
    # source unknown (this isn't used anywhere!)
    SRC_UKN = 0x10


# DQ FLAGS
class DQFlags(int, Enum):
    DQ_DO_NOT_USE = 0x01
    DQ_SATURATED = 0x02
    DQ_JUMP_DET = 0x04


# e name common names
SCI: Final[str] = "SCI"
BGD: Final[str] = "BGD"
RES: Final[str] = "RES"


# test states
class ExitStates(int, Enum):
    EXIT_SUCCESS = 0
    EXIT_FAIL = 1
    EXIT_EARLY = 2
    EXIT_MIXED = 3


# table column enum to be used to match table col names
# noinspection SpellCheckingInspection
class TableColumn(str, Enum):
    """Table column names used across the pipeline."""
    CAT_NUM = "Catalogue_Number"
    RA = "RA"
    DEC = "DEC"
    FLUX = "flux"
    E_FLUX = "eflux"
    FLUX_2 = "flux_2"
    X_CENTROID = "x_centroid"
    Y_CENTROID = "y_centroid"
    X_PEAK = "x_peak"
    Y_PEAK = "y_peak"
    EE_FRACTION = "eefraction"
    RADIUS = "radius"
    AP_CORR = "apcorr"
    STD_FLUX = "stdflux"
    NUM = "NUM"
    FLAG = "flag"
    FLUX_DET = "flux_det"
    FLUX_FIT = "flux_fit"
    FLUX_ERR = "flux_err"
    OUT_FLUX = "outflux"
    X_0 = "x_0"
    Y_0 = "y_0"
    X_DET = "x_det"
    Y_DET = "y_det"
    ID = "id"
    MAG = "mag"
    MAG_UPPER = "MAG"
    MAG_DET = "mag_det"
    MAG_DIFF = "mag_diff"
    ERROR_MAG = "eMAG"
    STATUS = "status"
    REC = "rec"
    PARAM = "PARAM"
    X_INIT = "x_init"
    Y_INIT = "y_init"
    XY_DEV = "xydev"
    XY_DEV_ = "_xydev"
    ERR_LOWER = "err"
    OFF = "off"
    X_FIT = "x_fit"
    Y_FIT = "y_fit"
    Q_FIT = "qfit"
    PUPIL = "pupil"
    SKY = "sky"
    SMOOTHNESS = "smoothness"
    SHARPNESS = "sharpness"
    ROUNDNESS1 = "roundness1"
    ROUNDNESS2 = "roundness2"
    RA_1 = "RA_1"
    RA_2 = "RA_2"
    TOTAL_FLUX_ADDED = "TOTAL_FLUX_ADDED"
    TOTAL_MAG = "TOTAL_MAG"
    TOTAL_DIFF_MAG = "TOTAL_DIFF_MAG"

    # needed as the table system doenst seem to handle enums properly
    def __str__(self) -> str:
        return self.value

    # needed as the table system doenst seem to handle enums properly
    def __format__(self, format_spec: str) -> str:
        return self.value.__format__(format_spec)


# DEFAULT MATCHING COLS
MATCH_COLS: List[str] = [
    TableColumn.RA, TableColumn.DEC, TableColumn.FLAG, TableColumn.FLUX,
    TableColumn.E_FLUX, TableColumn.NUM]


# Q table col names
class QTableColNames(str, Enum):
    SUM_ERR_0 = "aperture_sum_err_0"
    SUM_0 = "aperture_sum_0"
    SUM_1 = "aperture_sum_1"

    # needed as the table system does not seem to handle enums properly
    def __str__(self) -> str:
        return self.value

    # needed as the table system does not seem to handle enums properly
    def __format__(self, format_spec: str) -> str:
        return self.value.__format__(format_spec)


# tag for header
# noinspection SpellCheckingInspection
class HeaderTags(str, Enum):
    FILTER_LOWER = "filter"
    FILTER = "FILTER"
    EXT = "XTENSION"
    IMAGE = "IMAGE"
    BIN_TABLE = "BINTABLE"
    OUTPUT = "OUTPUT"
    STAR_BUG = "STARBUG"
    CALIBRATION_LV = "CALIBLEVEL"
    NAXIS = "NAXIS"
    NAXIS1 = "NAXIS1"
    NAXIS2 = "NAXIS2"
    C_TYPE = "CTYPE"
    OBS = "OBSERVTN"
    VISIT = "VISIT"
    EXPOSURE = "EXPOSURE"

    # needed as the table system doenst seem to handle enums properly
    def __str__(self) -> str:
        return self.value

    # needed as the table system doenst seem to handle enums properly
    def __format__(self, format_spec: str) -> str:
        return self.value.__format__(format_spec)


# tags for image header
# noinspection SpellCheckingInspection
class ImageHeaderTags(str, Enum):
    DETECTOR = "DETECTOR"
    TELESCOPE = "TELESCOP"
    INSTRUMENT = "INSTRUME"
    BUN_IT = "BUNIT"
    PIXAR_A2 = "PIXAR_A2"
    PIXAR_SR = "PIXAR_SR"
    JWST = "JWST"
    FILTER = "FILTER"

    # needed as the table system doenst seem to handle enums properly
    def __str__(self) -> str:
        return self.value

    # needed as the table system doenst seem to handle enums properly
    def __format__(self, format_spec: str) -> str:
        return self.value.__format__(format_spec)


# tag used for param file.
VERBOSE_TAG: Final[str] = "VERBOSE"


# mode labels.
# noinspection SpellCheckingInspection
class Modes(str, Enum):
    DETECTION = "DETECTION"
    BACKGROUND = "BACKGROUND"
    APP_HOT = "APPHOT"
    PSFP_HOT = "PSFPHOT"
    MATCH_OUTPUTS = "MATCHOUTPUTS"
    CLEAR = "CLEAR"


# HASH DEFS
STAR_BUG_MIRI: Final[int] = 1
NIRCAM: Final[int] = 2
NIRCAM_STRING: Final[str] = "NIRCAM"
MIRI_STRING: Final[str] = "MIRI"
# noinspection SpellCheckingInspection
MIRI_IMAGE = "MIRIMAGE"


class DetectorLengths(int, Enum):
    NULL = 0
    LONG = 1
    SHORT = 2


# enum unit
# noinspection SpellCheckingInspection
class Units(int, Enum):
    PIX = 0
    ARCSEC = 1
    ARCMIN = 2
    DEG = 3


# the column names of the rsult table used by artificial stars
TEST_TABLE_COLUMN_NAMES: Final[List[str]] = [
    TableColumn.X_0, TableColumn.Y_0, TableColumn.MAG, TableColumn.FLUX,
    TableColumn.X_DET, TableColumn.Y_DET, TableColumn.FLUX_DET,
    TableColumn.STATUS]


# text based logo (using raw string to bypass escape characters)
_LOGO_PATH = Path(os.path.join(
    os.path.join(Path(__file__).parent, "extras"), "logo.txt"))
LOGO: Final[str] = _LOGO_PATH.read_text(encoding="utf-8") + "%s"

# dictionary of help strings for specific modes (
# DETECTION, BACKGROUND, APP_HOT, PSF_PHOT, MATCH_OUTPUTS).
# noinspection SpellCheckingInspection
HELP_STRINGS = {
    Modes.DETECTION:
        """
            Source Detection
            ----------------

            This routine locates point sources in an image. The input is a
            FITS image and the output is a FITS table, containing a list of
            point source locations, their geometric properties and
            flux/magnitude measurements as calculated by aperture photometry.
            The output file will have the suffix "-ap", note this is the same
            as the output for the aperture photometry routine.

            To run this routine, use the core command:

                $~ starbug2 -D image.fits

            Alter the parameter file options under "DETECTION" to tune the
            performance of starbug2. Two of the key parameters are:

                - SIGSKY : Set the background level of the image
                - SIGSRC : Set the detection threshold of the sources

            Full documentation is at https://starbug2.readthedocs.io
        """,
    Modes.BACKGROUND:
        """
            Diffuse Background Estimations
            ------------------------------

            This routine estimates the "dusty" emissions in an image, given
            a source list. It is used to subtract from the image, thus removing
            the flux contribution on a source brightness from the dusty
            environment.

            The routine requires a list of sources to be generated (by source
            detection) or loaded with [-d sourcelist.fits] and requires a FITS
            image to work on. The routine will ouput a FITS image, with the
            same dimensions and spatial coverage as the input image, with the
            suffix "-bgd". This background image can be used in the photometry
             later.

            To run the routine, use the core command:

                $~ starbug2 -B -d sourcelist.fits image.fits

            Alter the parameter file options under "BACKGROUND ESTIMATION"
            to tune the performance of starbug2. Two key parameters are:

                - BGD_R    : Set a fixed aperture mask radius around each
                             source
                - BOX_SIZE : Set the estimation resolution (larger will be
                             more blurred)

            Full documentation is at https://starbug2.readthedocs.io
        """,
    Modes.APP_HOT:
        """
            Aperture Photometry
            -------------------

            This routine conducts aperture photometry on an image given a list
            of sources. It requires a FITS image to run on and a FITS table
            source list with either RA/DEC columns, or x/y_centroid or x/y_0
            columns. The routine outputs a table with the suffix "-ap". Note
            this filename is the same as the source detection routine because
            aperture photometry is automatically run at the end of the source
            detection step. The output table contains 2flux/magnitude
            information on every source

            To run this routine, use the core command:

                $~ starbug2 -A -d sourcelist.fits image.fits

            Alter the parameter file options under "APERTURE PHOTOMETRY" to
            tune the performance of starbug2. Three key parameters are:

                - APPHOT_R : Set the aperture radius for photometry (in pixels)
                - SKY_RIN  : Set the inner sky annulus radius (in pixels)
                - SKY_ROUT : Set the outer sky annulus radius (in pixels)

            Full documentation is at https://starbug2.readthedocs.io
        """,
    Modes.PSFP_HOT:
        """
            PSF Photometry
            --------------

            This routine conducts PSF fitting photometry on an image given
            a list of sources. Its requires a FITS image to run on and a FITS
            table sourcelist with either RA/DEC columns, or x/y_centroid or
            x/y_0 columns. The routine outputs a table with the suffix "-psf".
            The output table contains 2flux/magnitude information on every
            source.

            To run this routine, use the core command:

                $~ starbug2 -P -d sourcelist image.fits

            Alter the parameter file options under "PHOTOMETRY" to tune the
            performance of starbug2. Two key parameters are:

                - FORCE_POS    : Hold the cetroid positions of source fixed
                                 (forced photometry)
                - GEN_RESIDUAL : Generate a residual image from all the fit
                                 source

            Full documentation is at https://starbug2.readthedocs.io
        """,
    Modes.MATCH_OUTPUTS:
        """
            Match Outputs
            -------------

            This option is set if the user wishes to combine all the output
            catalogues from starbug together. It would be used in the case
            that a routine is being ran on a list of images (either in series
            or parallel) and the final catalogues should all be combined into
            a single source list. It outputs two files, one with the suffix
            "full" and another with "match". The first is all columns from all
            table preserved into a single large catalogue, the second averages
            all the similar columns into a reduced table.

            To run this routine, use the core code:

                $~ starbug2 -DM image1.fits image2.fits image3.fits ...

            Alter the parameter file options under "CATALOGUE MATCHING" to
            tune the performance of starbug2. Two key parameters are:

                - MATCH_THRESH : Set the separation threshold (arcsec) to match
                                 two sources
                - NEXP_THRESH  : Set the minimum number of catalogues a source
                                 must be present in
        """,
}
