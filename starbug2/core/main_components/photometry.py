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
from typing import List, Callable, Tuple, Any, cast

import numpy as np
from astropy.io.fits import (
    HDUList, open, Header, ImageHDU, BinTableHDU, PrimaryHDU)
from astropy.stats import sigma_clipped_stats
from astropy.table import Table, Column, vstack
from astropy.wcs import WCS
from photutils.datasets import make_model_image

from photutils.psf import ImagePSF
from starbug2.constants import (
    TableColumn, ExitStates, HeaderTags, SourceFlags, Units, ImageHeaderTags,
    MIRI_STRING, MIRI_IMAGE, DetectorLengths, NIRCAM, BGD_FILE, AP_FILE)
from starbug2.core.star_bug_config import StarBugMainConfig
from starbug2.interfaces.star_bug_interface import StarBugInterface
from starbug2.routines.psf_phot_routine import PSFPhotRoutine
from starbug2.utilities.filters import STAR_BUG_FILTERS, FilterStruct
from starbug2.utilities.utils import (
    flux_to_pogson_mag, p_error, parse_unit, warn, get_data_path,
    get_mj_ysr2jy_scale_factor, reindex)


class Photometry:

    @staticmethod
    def _determine_f_name(filter_string: str | None, f_name: str | None,
                          info: dict[str, str]) -> str | None:
        """
        Finds the file name for the PSF if not provided
        :param filter_string: Name of the photometric filter band used.
        :type filter_string: str | None
        :param f_name: Filename of a PSF fits image
        :type f_name: str
        :param info: Metadata dictionary detailing execution properties.
        :type info: dict[str, str]
        :return: the correct psf file name
        :rtype: str
        """
        assert filter_string is not None
        if not f_name:
            filter_struct: FilterStruct | None = (
                STAR_BUG_FILTERS.get(filter_string))
            if filter_struct:
                dt_name: str = info[ImageHeaderTags.DETECTOR]
                # noinspection SpellCheckingInspection
                if dt_name == "NRCALONG":
                    # noinspection SpellCheckingInspection
                    dt_name = "NRCA5"
                # noinspection SpellCheckingInspection
                if dt_name == "NRCBLONG":
                    # noinspection SpellCheckingInspection
                    dt_name = "NRCB5"
                if dt_name == "MULTIPLE":
                    if (filter_struct.instr == NIRCAM
                            and filter_struct.length == DetectorLengths.SHORT):
                        # noinspection SpellCheckingInspection
                        dt_name = "NRCA1"
                    elif (filter_struct.instr == NIRCAM and
                          filter_struct.length == DetectorLengths.LONG):
                        # noinspection SpellCheckingInspection
                        dt_name = "NRCA5"
                    elif filter_struct.instr == MIRI_STRING:
                        dt_name = ""
                if dt_name == MIRI_IMAGE:
                    dt_name = ""
                f_name = "%s/%s%s.fits" % (
                    get_data_path(), filter_string, dt_name)
            else:
                f_name = None
        return f_name

    # noinspection SpellCheckingInspection
    @staticmethod
    def load_psf(
            filter_string: str | None, info: dict[str, str],
            log: Callable[[str], None],
            f_name: str | None = None) -> Tuple[ExitStates, np.ndarray | None]:
        """
        Load a PSF_FILE to be used during photometry.
        :param filter_string: Name of the photometric filter band used.
        :type filter_string: str | None
        :param info: Metadata dictionary detailing execution properties.
        :type info: dict[str, str]
        :param log: Callable logging function for status tracking.
        :type log: Callable[[str], None]
        :param f_name: Filename of a PSF fits image
        :type f_name: str
        :return: The exit status state and the psf
        :rtype: Tuple[ExitStates, np.ndarray | None]
        """
        status: ExitStates = ExitStates.EXIT_SUCCESS
        f_name = Photometry._determine_f_name(filter_string, f_name, info)

        psf: np.ndarray | None = None
        if f_name is not None and os.path.exists(f_name):
            fp: HDUList = open(f_name)

            if fp[0].data is None:
                (status, psf) = Photometry.load_psf(
                    filter_string, info, log, None)
                if status != ExitStates.EXIT_SUCCESS:
                    p_error(
                        "There is a version mismatch between starbug and "
                        "webbpsf. Please reinitialise with: "
                        "starbug2 --init.\n")
                    quit("Fatal error, quitting\n")
                return status, psf
            psf = fp[0].data
            fp.close()
            log("loaded PSF_FILE='%s'\n" % f_name)
        else:
            p_error("PSF_FILE='%s' does not exist\n" % f_name)
            status = ExitStates.EXIT_FAIL
        return status, psf

    @staticmethod
    def input_checks(
            filter_string: str | None,
            main_image: ImageHDU | PrimaryHDU | None,
            detections: Table | None,
            psf: np.ndarray | None,
            wcs: WCS | None, info, log, config) -> (
                Tuple[ExitStates, np.ndarray | None]):
        """
        Validates pipeline input structures prior to running photometry.

        Verifies the presence of mandatory structures including the image,
        coordinate system, and source catalogue. Automatically attempts to
        resolve and load a valid PSF model if one is not already provided.

        :param filter_string: Name of the photometric filter band used.
        :type filter_string: str | None
        :param main_image: Primary target image array used for measurements.
        :type main_image: ImageHDU | PrimaryHDU | None
        :param detections: Table of source positions or None if empty.
        :type detections: Table | None
        :param psf: Empirical or modelled point spread function array.
        :type psf: np.ndarray | None
        :param wcs: World Coordinate System object for spatial alignment.
        :type wcs: WCS | None
        :param info: Metadata dictionary detailing execution properties.
        :type info: dict[str, str]
        :param log: Callable logging function for status tracking.
        :type log: Callable[[str], None]
        :param config: Configuration instance containing engine parameters.
        :type config: StarBugMainConfig
        :return: A tuple containing the execution exit status state and the
                 resolved PSF model array (or None if validation failed).
        :rtype: Tuple[ExitStates, np.ndarray | None]
        """
        if filter_string is None or wcs is None:
            return ExitStates.EXIT_FAIL, None
        if main_image is None:
            return ExitStates.EXIT_FAIL, None
        if detections is None:
            p_error("unable to run photometry: no source list loaded\n")
            return ExitStates.EXIT_FAIL, None

        if psf is None:
            result: ExitStates
            (result, psf) = Photometry.load_psf(
                filter_string, info, log,
                os.path.expandvars(config.psf_file_override))
            if result != ExitStates.EXIT_SUCCESS:
                p_error("unable to run photometry: no PSF loaded\n")
                return ExitStates.EXIT_FAIL, None
        return ExitStates.EXIT_SUCCESS, psf

    @staticmethod
    def prepare_data(
            image: HDUList | None, log: Callable[[str], None],
            background: ImageHDU | PrimaryHDU | None, header: Header,
            main_image: ImageHDU | PrimaryHDU,
            config: StarBugMainConfig, full_width_half_max: float) -> (
                Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]):
        """
        Prepares and conditions image arrays for pipeline processing.

        Initialises primary data, error, and mask arrays. If no external
        background map is provided, a robust background level is measured
        on-the-fly using sigma-clipped statistics.

        :param image: Full FITS HDU list containing the data extensions.
        :type image: HDUList | None
        :param log: Callable logging function for status tracking.
        :type log: Callable[[str], None]
        :param background: Estimated background map array image if loaded.
        :type background: ImageHDU | PrimaryHDU | None
        :param header: FITS header containing primary data metadata.
        :type header: Header
        :param main_image: Primary target image array used for measurements.
        :type main_image: ImageHDU | PrimaryHDU
        :param config: Configuration instance containing engine parameters.
        :type config: StarBugMainConfig
        :param full_width_half_max: the full width 1/2 max value.
        :type full_width_half_max: float
        :return: A tuple containing four fully processed NumPy arrays:
                 (image_data, error, background, mask, min_separation).
        :rtype: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]
        """
        image_data: np.ndarray
        error: np.ndarray
        bgd: np.ndarray | None
        mask: np.ndarray
        image_data, error, bgd, mask = (
            StarBugInterface.prepare_image_arrays(
                image, log, background, header, main_image))

        if bgd is None:
            clipped_median: float
            _, clipped_median, _ = (
                sigma_clipped_stats(image_data, sigma=config.sigma_sky))
            bgd = cast(np.ndarray, np.ones(main_image.shape) * clipped_median)
            log(
                "-> no background file loaded, measuring sigma "
                "clipped median\n")
        assert bgd is not None

        min_separation: float = config.critical_separation
        if not min_separation:
            min_separation = min(5.0, 2.5 * full_width_half_max)

        return image_data, error, bgd, mask, min_separation

    @staticmethod
    def create_psf_model_mask(
            psf: np.ndarray | None, config: StarBugMainConfig,
            log: Callable[[str], None]) -> Tuple[np.ndarray, ImagePSF, int]:
        """
        Cleans the PSF array, wraps it in a model, and calculates fit size.

        Identifies and masks out non-finite pixels (NaN/Inf) within the raw
        PSF data, builds an ImagePSF instance, and enforces an odd integer
        pixel width constraint for the profiling boundary.

        :param psf: Raw empirical or modelled point spread function array.
        :type psf: np.ndarray
        :param config: Configuration instance containing engine parameters.
        :type config: StarBugMainConfig
        :param log: Callable logging function for status tracking.
        :type log: Callable[[str], None]
        :return: A tuple containing the updated boolean mask array, the
                 initialized ImagePSF model object, and the computed odd-valued
                 fitting bounding dimension size.
        :rtype: Tuple[np.ndarray, ImagePSF, int]
        """
        psf_mask: np.ndarray = ~np.isfinite(psf)
        if psf_mask.sum():
            assert psf is not None
            psf[psf_mask] = 0
            log("-> masking INF pixels in PSF_FILE\n")

        psf_model: ImagePSF = ImagePSF(data=psf)
        size: int
        if config.psf_fit_size > 0:
            size = config.psf_fit_size
        else:
            size = psf_model.shape[0]
        if not size % 2:
            size -= 1
        log("-> psf size: %d\n" % size)
        return psf_mask, psf_model, size

    @staticmethod
    def determine_initial_guesses(
            config: StarBugMainConfig, detections: Table,
            main_image: ImageHDU | PrimaryHDU, filter_string: str) -> (
                Tuple[float, Table]):
        """
        Extracts starting positions and sets the target aperture radius.

        Parses the configuration parameters to isolate the designated initial
        aperture profiling width and extracts source detection columns
        associated with coordinate positions and filter-specific fluxes to
        seed the optimisation engine.

        :param config: Configuration instance containing engine parameters.
        :type config: StarBugMainConfig
        :param detections: Table containing detected sources and properties.
        :type detections: Table
        :param main_image: Primary target image extension containing data
                           units.
        :type main_image: ImageHDU | PrimaryHDU
        :param filter_string: Name of the photometric filter band used.
        :type filter_string: str
        :return: A tuple containing the evaluated aperture hot-radius value
                 and an astropy Table tracking input parameter estimates.
        :rtype: Tuple[float, Table]
        """

        # Sort out Init guesses
        app_hot_r: float = float(config.aperture_phot_radius)
        if not app_hot_r or app_hot_r <= 0:
            app_hot_r = 3.0

        init_guesses: Table = detections.copy()

        # noinspection DuplicatedCode
        if TableColumn.X_CENTROID in init_guesses.colnames:
            init_guesses.rename_column(
                TableColumn.X_CENTROID, TableColumn.X_INIT)
        if TableColumn.Y_CENTROID in init_guesses.colnames:
            init_guesses.rename_column(
                TableColumn.Y_CENTROID, TableColumn.Y_INIT)
        if TableColumn.X_DET in init_guesses.colnames:
            init_guesses.rename_column(
                TableColumn.X_DET, TableColumn.X_INIT)
        if TableColumn.Y_DET in init_guesses.colnames:
            init_guesses.rename_column(
                TableColumn.Y_DET, TableColumn.Y_INIT)

        init_guesses = init_guesses[init_guesses[TableColumn.X_INIT] >= 0]
        init_guesses = init_guesses[init_guesses[TableColumn.Y_INIT] >= 0]
        init_guesses = init_guesses[
            init_guesses[TableColumn.X_INIT]
            < main_image.header[HeaderTags.NAXIS1]]
        init_guesses = init_guesses[
            init_guesses[TableColumn.Y_INIT]
            < main_image.header[HeaderTags.NAXIS2]]

        # Allow tables that don't have the correct columns through
        required: List[str] = [
            TableColumn.X_INIT, TableColumn.Y_INIT, TableColumn.FLUX,
            filter_string, TableColumn.FLAG]
        for notfound in set(required) - set(init_guesses.colnames):
            dtype = np.uint16 if notfound == TableColumn.FLAG else float
            init_guesses.add_column(
                Column(np.zeros(len(init_guesses)),
                       name=notfound, dtype=dtype))

        init_guesses = init_guesses[required]
        init_guesses.remove_column(TableColumn.FLUX)
        init_guesses.rename_column(filter_string, "ap_%s" % filter_string)

        return app_hot_r, init_guesses

    @staticmethod
    def _convert_to_pixel(header: Header, max_y_dev: float) -> float:
        """
        Converts an angular value in arc seconds to its pixel equivalent.

        Extracts the pixel area scale factor from the image header to evaluate
        the linear pixel scale, then divides the input arcsecond displacement
        value to project it into pixel units. Logs a warning if the scale
        metadata is missing.

        :param header: FITS header metadata housing the projection parameters.
        :type header: Header
        :param max_y_dev: Angular maximum displacement threshold in arc
                          seconds.
        :type max_y_dev: float
        :return: The resolved tracking displacement limit expressed in pixels.
        :rtype: float
        """
        if not header.get(ImageHeaderTags.PIXAR_A2):
            # noinspection SpellCheckingInspection
            warn(
                "MAX_XYDEV is units arcseconds, but starbug "
                "cannot locate a pixel scale in the header."
                " Please use syntax MAX_XYDEV=%sp to set "
                "change to pixels\n" % max_y_dev)
        else:
            max_y_dev /= np.sqrt(
                header.get(ImageHeaderTags.PIXAR_A2))
        return max_y_dev

    @staticmethod
    def convert_unit_to_max_y_dev(
            max_y_dev: float, unit: int, header: Header) -> float:
        """
        Converts a maximum deviation threshold value to pixel units.

        Parses astronomical tracking limits expressed across varying spatial
        units (Degrees, Arc minutes, Arc seconds, or raw Pixels), scales the
        amplitudes appropriately, and resolves angular values to pixel
        dimensions using the image header plate scale factors.

        :param max_y_dev: Input maximum displacement tracking limit value.
        :type max_y_dev: float
        :param unit: Enumerated identifier representing the measurement scale.
        :type unit: int
        :param header: FITS header metadata housing the projection parameters.
        :type header: Header
        :return: Computed spatial maximum threshold value expressed in pixels.
        :rtype: float
        """
        if unit is not None:
            match unit:
                case Units.DEG:
                    max_y_dev *= 3600
                    max_y_dev = Photometry._convert_to_pixel(header, max_y_dev)
                case Units.ARCMIN:
                    max_y_dev *= 60
                    max_y_dev = Photometry._convert_to_pixel(header, max_y_dev)
                case Units.ARCSEC:
                    max_y_dev = Photometry._convert_to_pixel(header, max_y_dev)
                case Units.PIX:
                    pass
                case _:
                    p_error("Dont recognise the unit.")
        return max_y_dev

    @staticmethod
    def execute_revalidation_of_positions(
            psf_cat: Table, image_data: np.ndarray, error: np.ndarray,
            mask: np.ndarray, psf_model: ImagePSF, size: int,
            min_separation: float, app_hot_r: float, bgd: np.ndarray,
            filter_string: str, config: StarBugMainConfig,
            log: Callable[[str], None], max_y_dev: float):
        """
        Executes forced-position profiling for astrometric outlier sources.

        Isolates individual sources whose fitted coordinates drift beyond
        the maximum permissible pixel variation threshold, re-evaluates them
        by locking their centroids to initial coordinates, updates diagnostic
        flags, and integrates the results back into the primary catalogue.

        :param psf_cat: Table containing calibrated stellar fit properties.
        :type psf_cat: Table
        :param image_data: Cleaned multidimensional science data array.
        :type image_data: np.ndarray
        :param error: Estimated variance uncertainties or weight maps.
        :type error: np.ndarray
        :param mask: Boolean flags highlighting pixels to skip.
        :type mask: np.ndarray
        :param psf_model: Initialised empirical or analytical PSF model.
        :type psf_model: ImagePSF
        :param size: Dimension bounding length of the fitting footprint box.
        :type size: int
        :param min_separation: Minimum pixel separation between distinct
                               sources.
        :type min_separation: float
        :param app_hot_r: Base aperture radius value evaluated for fits.
        :type app_hot_r: float
        :param bgd: Evaluated background grid or uniform pixel array.
        :type bgd: np.ndarray
        :param filter_string: Name of the photometric filter band used.
        :type filter_string: str
        :param config: Configuration instance containing engine parameters.
        :type config: StarBugMainConfig
        :param log: Callable status monitoring and message log tracker.
        :type log: Callable[[str], None]
        :param max_y_dev: Maximum allowable coordinate deviation threshold in
                          pixels.
        :type max_y_dev: float
        :return: Updated source catalogue containing the fixed outlier
                 measurements.
        :rtype: Table
        """
        if max_y_dev > 0:
            log(
                "-> position fit threshold: %.2gpix\n" % max_y_dev)
            phot = PSFPhotRoutine(
                psf_model, size, min_separation=min_separation,
                app_hot_r=app_hot_r, background=bgd, force_fit=1,
                verbose=config.verbose_logs)
            ii: np.ndarray = psf_cat[TableColumn.XY_DEV] > max_y_dev
            fixed_centres: Table = psf_cat[ii][
                [TableColumn.X_INIT, TableColumn.Y_INIT,
                 "ap_%s" % filter_string, TableColumn.FLAG]]
            if len(fixed_centres):
                log("-> forcing positions for deviant sources\n")
                fixed_cat: Table = phot(
                    image_data, init_params=fixed_centres,
                    error=error, mask=mask)
                # ABS. why are we using such an aggressive type check
                # here?
                fixed_cat[TableColumn.FLAG] |= (
                    np.uint16(SourceFlags.SRC_FIX))
                psf_cat.remove_rows(ii.tolist())
                psf_cat = vstack((psf_cat, fixed_cat))
            else:
                log("-> no deviant sources\n")
        return psf_cat

    @staticmethod
    def revalidate_deviant_positions(
            psf_cat: Table, image_data: np.ndarray, error: np.ndarray,
            mask: np.ndarray, psf_model: ImagePSF, size: int,
            min_separation: float, app_hot_r: float, bgd: np.ndarray,
            filter_string: str, header: Header,
            config: StarBugMainConfig, log: Callable[[str], None]) -> Table:
        """
        Identifies and re-fits sources exceeding centroid deviation thresholds.

        Parses spatial deviation tolerances across distinct coordinate units,
        isolates outlier astrometric fits via comparison to original tracking
        positions, and enforces a static position constraint profile fit onto
        anomalous targets to prevent mathematical divergence.

        :param psf_cat: Table containing calibrated stellar fit properties.
        :type psf_cat: Table
        :param image_data: Cleaned multidimensional science data array.
        :type image_data: np.ndarray
        :param error: Estimated variance uncertainties or weight maps.
        :type error: np.ndarray
        :param mask: Boolean flags highlighting pixels to skip.
        :type mask: np.ndarray
        :param psf_model: Initialised empirical or analytical PSF model.
        :type psf_model: ImagePSF
        :param size: Dimension bounding length of the fitting footprint box.
        :type size: int
        :param min_separation: Minimum pixel separation between distinct
                               sources.
        :type min_separation: float
        :param app_hot_r: Base aperture radius value evaluated for fits.
        :type app_hot_r: float
        :param bgd: Evaluated background grid or uniform pixel array.
        :type bgd: np.ndarray
        :param filter_string: Name of the photometric filter band used.
        :type filter_string: str
        :param header: Structural FITS header metadata for mapping updates.
        :type header: Header
        :param config: Configuration instance containing engine parameters.
        :type config: StarBugMainConfig
        :param log: Callable status monitoring and message log tracker.
        :type log: Callable[[str], None]
        :return: Updated source catalogue tracking resolved stellar parameters.
        :rtype: Table
        """
        # Setting position max variation
        max_y_dev: float
        unit: int
        max_y_dev, unit = parse_unit(config.max_xy_deviation)

        if unit is None or max_y_dev <= 0:
            return psf_cat

        max_y_dev = Photometry.convert_unit_to_max_y_dev(
            max_y_dev, unit, header)
        psf_cat = Photometry.execute_revalidation_of_positions(
            psf_cat, image_data, error, mask, psf_model, size, min_separation,
            app_hot_r, bgd, filter_string, config, log, max_y_dev)
        return psf_cat

    @staticmethod
    def export_psf_catalogue(
            psf_cat: Table, wcs: WCS, filter_string: str | None,
            header: Header, config: StarBugMainConfig, ap_file: str | None,
            background_file: str | None, out_dir: str | None,
            b_name: str | None, log: Callable[[str], None]) -> Table:
        """
        Calculates physical properties, appends metadata, and saves table.

        Transforms raw pixel positions to celestial coordinates ($RA, Dec$)
        using the WCS object, converts fitted instrumental fluxes to absolute
        magnitudes, aggregates diagnostic file paths into table metadata,
        and writes out a finalised FITS binary table.

        :param psf_cat: Table containing calibrated stellar fit properties.
        :type psf_cat: Table
        :param wcs: World Coordinate System conversion instance.
        :type wcs: WCS
        :param filter_string: Name of the photometric filter band used.
        :type filter_string: str | None
        :param header: Structural FITS header metadata for mapping updates.
        :type header: Header
        :param config: Configuration instance containing engine parameters.
        :type config: StarBugMainConfig
        :param ap_file: File path to the aperture correction table used.
        :type ap_file: str | None
        :param background_file: File path to the background map used.
        :type background_file: str | None
        :param out_dir: Target file-system directory paths for exported tables.
        :type out_dir: str | None
        :param b_name: Base filename identifier prefix for outputs.
        :type b_name: str | None
        :param log: Callable status monitoring and message log tracker.
        :type log: Callable[[str], None]
        :return: The finalised, fully populated astronomical source catalogue.
        :rtype: Table
        """
        ra: np.ndarray
        dec: np.ndarray
        ra, dec = wcs.all_pix2world(
            psf_cat[TableColumn.X_FIT], psf_cat[TableColumn.Y_FIT], 0)
        psf_cat.add_column(Column(ra, name=TableColumn.RA), index=2)
        psf_cat.add_column(Column(dec, name=TableColumn.DEC), index=3)

        mag: float
        mag_err: float
        mag, mag_err = flux_to_pogson_mag(
            psf_cat[TableColumn.FLUX], psf_cat[TableColumn.E_FLUX])

        filter_string: str = (
            filter_string if filter_string else TableColumn.MAG)
        psf_cat.add_column(
            mag + config.zero_point_magnitude, name=filter_string)
        psf_cat.add_column(mag_err, name="e%s" % filter_string)
        psf_catalogue = psf_cat

        # verify catalogue isn't none
        assert psf_catalogue is not None

        psf_catalogue.meta = dict(header.items())
        psf_catalogue.meta[AP_FILE] = ap_file
        psf_catalogue.meta[BGD_FILE] = background_file

        reindex(psf_catalogue)

        file_name: str = (
            "%s/%s-psf.fits" % (out_dir, b_name))
        log("--> %s\n" % file_name)
        BinTableHDU(
            data=psf_catalogue,
            header=header).writeto(file_name, overwrite=True)
        return psf_cat

    @staticmethod
    def generate_residual_image(
            psf_cat: Table, image_data: np.ndarray, psf_model: ImagePSF,
            size: int, bgd: np.ndarray, main_image: ImageHDU | PrimaryHDU,
            header: Header, wcs: WCS, out_dir: str | None,
            b_name: str | None, config: StarBugMainConfig,
            log: Callable[[str], None]) -> np.ndarray | None:
        # noinspection SpellCheckingInspection
        """
        Creates a source-subtracted residual map and writes it to disk.

        Subtracts the evaluated background and the fitted point spread function
        (PSF) stellar models from the target science data. Converts the image
        units to Janskys via scaling factors, injects updated WCS structural
        metadata into the header, and exports the map to a FITS file.

        :param psf_cat: Table containing calibrated stellar fit properties.
        :type psf_cat: Table
        :param image_data: Cleaned multidimensional science data array.
        :type image_data: np.ndarray
        :param psf_model: Initialised empirical or analytical PSF model.
        :type psf_model: ImagePSF
        :param size: Dimension bounding length of the fitting footprint box.
        :type size: int
        :param bgd: Evaluated background grid or uniform pixel array.
        :type bgd: np.ndarray
        :param main_image: Primary target image extension containing unit
                           metadata.
        :type main_image: ImageHDU | PrimaryHDU
        :param header: Structural FITS header metadata for mapping updates.
        :type header: Header
        :param wcs: World Coordinate System conversion instance.
        :type wcs: WCS
        :param out_dir: Target file-system directory paths for exported maps.
        :type out_dir: str | None
        :param b_name: Base filename identifier prefix for outputs.
        :type b_name: str | None
        :param config: Configuration instance containing engine parameters.
        :type config: StarBugMainConfig
        :param log: Callable status monitoring and message log tracker.
        :type log: Callable[[str], None]
        :return: Unit-scaled residual image array, or None if disabled.
        :rtype: np.ndarray | None
        """
        residuals: np.ndarray | None = None
        if config.generate_residual_image:
            log("-> generating residual\n")
            _tmp: Table = psf_cat[
                TableColumn.X_FIT, TableColumn.Y_FIT,
                TableColumn.FLUX].copy()
            _tmp.rename_columns(
                (TableColumn.X_FIT, TableColumn.Y_FIT),
                (TableColumn.X_0, TableColumn.Y_0))
            stars: np.ndarray = make_model_image(
                image_data.shape, psf_model, _tmp,
                model_shape=(size, size))

            residual: np.ndarray = image_data - (bgd + stars)
            residuals = (
                residual / get_mj_ysr2jy_scale_factor(main_image))
            header: Header = header
            header.update(wcs.to_header())
            ImageHDU(
                data=cast(Any, residuals),
                name="RES", header=header).writeto(
                "%s/%s-res.fits" % (out_dir, b_name),
                overwrite=True)
        return residuals

    @staticmethod
    def photometry_routine(
            filter_string: str | None, wcs: WCS | None,
            config: StarBugMainConfig,
            main_image: ImageHDU | PrimaryHDU | None,
            log: Callable[[str], None],
            image: HDUList | None,
            info: dict[str, str],
            background: ImageHDU | PrimaryHDU | None,
            header: Header,
            detections: Table | None,
            psf: np.ndarray | None,
            full_width_half_max: float,
            ap_file: str | None,
            background_file: str | None,
            out_dir: str | None,
            b_name: str | None) -> Tuple[int, Table | None, np.ndarray | None]:
        """
        Executes the master PSF-fitting photometry processing pipeline.

        Validates physical inputs, extracts cropped science matrices,
        initialises  background evaluation fields, builds or maps empirical
        stellar profile models, resolves positional deviations across
        coordinate boundaries, and updates structural header maps before
        exporting the finalised catalogues and residual maps to disk.

        :param filter_string: Name of the photometric filter band used.
        :type filter_string: str | None
        :param wcs: World Coordinate System conversion instance.
        :type wcs: WCS | None
        :param config: Configuration instance containing engine parameters.
        :type config: StarBugMainConfig
        :param main_image: Primary target image extension containing data
                           units.
        :type main_image: ImageHDU | PrimaryHDU | None
        :param log: Callable status monitoring and message log tracker.
        :type log: Callable[[str], None]
        :param image: Complete FITS structure housing science extensions.
        :type image: HDUList | None
        :param info: Dictionary containing target object and runtime metadata.
        :type info: dict[str, str]
        :param background: Explicit background structure or image extension.
        :type background: ImageHDU | PrimaryHDU | None
        :param header: Structural FITS header metadata for mapping updates.
        :type header: Header
        :param detections: Table containing detected sources and properties.
        :type detections: Table | None
        :param psf: Array containing the input point spread function template.
        :type psf: np.ndarray | None
        :param full_width_half_max: Full width at half maximum value of stars.
        :type full_width_half_max: float
        :param ap_file: File path to the aperture correction table used.
        :type ap_file: str | None
        :param background_file: File path to the background map used.
        :type background_file: str | None
        :param out_dir: Target file-system directory paths for exported
                        products.
        :type out_dir: str | None
        :param b_name: Base filename identifier prefix for outputs.
        :type b_name: str | None
        :return: Success state flag, the completed stellar source catalogue,
                 and the background-subtracted residual image array.
        :rtype: Tuple[int, Table | None, np.ndarray | None]
        """
        results: ExitStates
        (results, psf) = Photometry.input_checks(
            filter_string, main_image, detections, psf, wcs, info, log, config)
        if results != ExitStates.EXIT_SUCCESS:
            return ExitStates.EXIT_FAIL, None, None

        assert main_image is not None
        assert detections is not None
        assert filter_string is not None
        assert wcs is not None

        log("\nRunning PSF Photometry\n")

        # lock the types.
        image_data: np.ndarray
        error: np.ndarray
        bgd: np.ndarray
        mask: np.ndarray
        image_data, error, bgd, mask, min_separation = (
            Photometry.prepare_data(
                image, log, background, header, main_image, config,
                full_width_half_max))

        # Collect relevant files and data
        psf_mask: np.ndarray
        psf_model: ImagePSF
        size: int
        psf_mask, psf_model, size = Photometry.create_psf_model_mask(
            psf, config, log)

        app_hot_r: float
        init_guesses: Table
        app_hot_r, init_guesses = Photometry.determine_initial_guesses(
            config, detections, main_image, filter_string)

        # Run Fit
        phot: PSFPhotRoutine = PSFPhotRoutine(
            psf_model, size, min_separation=min_separation,
            app_hot_r=app_hot_r, background=bgd,
            force_fit=int(config.force_centroid_position),
            verbose=config.verbose_logs)
        psf_cat: Table = phot(
            image_data, init_params=init_guesses, error=error,
            mask=mask)

        # handle deviant positions if not in centroid positions.
        if config.force_centroid_position:
            psf_cat[TableColumn.FLAG] |= SourceFlags.SRC_FIX
        else:
            if not psf_cat:
                return ExitStates.EXIT_FAIL, None, None
            psf_cat = Photometry.revalidate_deviant_positions(
                psf_cat, image_data, error, mask, psf_model, size,
                min_separation, app_hot_r, bgd, filter_string, header,
                config, log)
        Photometry.export_psf_catalogue(
            psf_cat, wcs, filter_string, header, config, ap_file,
            background_file, out_dir, b_name, log)

        # Residual Image
        residuals: np.ndarray | None = Photometry.generate_residual_image(
            psf_cat, image_data, psf_model, size, bgd, main_image, header, wcs,
            out_dir, b_name, config, log)
        return ExitStates.EXIT_SUCCESS, psf_cat, residuals
