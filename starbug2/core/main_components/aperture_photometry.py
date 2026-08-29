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
from typing import Callable, Tuple

import numpy as np
from astropy.io.fits import ImageHDU, PrimaryHDU, Header, HDUList
from astropy.table import Column, Table, hstack
from starbug2.constants import (
    ExitStates, TableColumn, ImageHeaderTags, NIRCAM_STRING, MIRI_STRING, DQ,
    HeaderTags)
from starbug2.core.star_bug_config import StarBugMainConfig
from starbug2.interfaces.star_bug_interface import StarBugInterface
from starbug2.routines.app_phot_routine import APPhotRoutine
from starbug2.utilities.utils import (
    p_error, warn, ext_names, flux_to_pogson_mag, reindex, export_table,
    get_data_path)


class AperturePhotometry:
    """
    this class contains the main logic for running aperture photometry from
    starbug main.
    """

    @staticmethod
    def check_states(detections, image, filter_string) -> ExitStates:
        """
        checks if the inputs are in a valid state.
        :param detections: Table of source positions or None if empty.
        :type detections: Table | None
        :param image: Full FITS HDU list containing the data extensions.
        :type image: HDUList | None
        :param filter_string: Name of the photometric filter band used.
        :type filter_string: str | None
        :return: exit fail if the inputs are not valid, success otherwise.
        :rtype: ExitStates
        """
        if detections is None:
            p_error("No detection source file loaded (-d file-ap.fits)\n")
            return ExitStates.EXIT_FAIL
        if image is None:
            p_error("No image provided")
            return ExitStates.EXIT_FAIL
        if len({TableColumn.X_0, TableColumn.Y_0, TableColumn.X_INIT,
                TableColumn.Y_INIT, TableColumn.X_CENTROID,
                TableColumn.Y_CENTROID} & set(detections.colnames)) < 2:
            p_error("No pixel coordinates in source file\n")
            return ExitStates.EXIT_FAIL
        if filter_string is None:
            p_error("no filter name")
            return ExitStates.EXIT_FAIL
        return ExitStates.EXIT_SUCCESS

    @staticmethod
    def remove_detection_columns(detections, filter_string) -> Table:
        """
        removes columns from the detections to prepare it for photometry.
        :param detections: Table of source positions
        :type detections: Table
        :param filter_string:  Name of the photometric filter band used.
        :type filter_string: str
        :return:
        """
        new_columns: tuple[str, str, str, str, str, str | None, str] = (
            TableColumn.SMOOTHNESS, TableColumn.FLUX, TableColumn.E_FLUX,
            TableColumn.SKY, TableColumn.FLAG, filter_string,
            "e%s" % filter_string)
        assert detections is not None
        detections.remove_columns(
            set(new_columns) & set(detections.colnames))
        return detections

    @staticmethod
    def aperture_correction(
            config: StarBugMainConfig, info: dict[str, str],
            log: Callable[[str], None], filter_string: str,
            full_width_half_max: float) -> Tuple[float, float, float, float]:
        """
        Determines the aperture radii and calculates the aperture correction.

        Resolves the appropriate instrument calibration file to compute the
        photometric aperture radius based on config constraints, falling back
        to profiling via FWHM or default values if required.

        :param config: Configuration instance containing engine parameters.
        :type config: StarBugMainConfig
        :param info: Metadata dictionary detailing execution properties.
        :type info: dict[str, str]
        :param log: Callable logging function for status tracking.
        :type log: Callable[[str], None]
        :param filter_string: Name of the photometric filter band used.
        :type filter_string: str
        :param full_width_half_max: FWHM value for stellar profiles.
        :type full_width_half_max: float
        :return: A tuple containing the calculated parameters: (aperture
                 radius, sky inner radius, sky outer radius, aperture
                 correction factor).
        :rtype: Tuple[float, float, float, float]
        """
        ap_corr_f_name: str | None = None
        if _ap_corr_f_name := config.ap_corr_file_override:
            ap_corr_f_name = _ap_corr_f_name
        elif info.get(ImageHeaderTags.INSTRUMENT) == NIRCAM_STRING:
            ap_corr_f_name = os.path.join(
                get_data_path(), "apcorr_nircam.fits")
        elif info.get(ImageHeaderTags.INSTRUMENT) == MIRI_STRING:
            ap_corr_f_name = os.path.join(
                get_data_path(), "apcorr_miri.fits")

        if ap_corr_f_name:
            log("-> apcorr file: %s\n" % ap_corr_f_name)
        else:
            warn("No apcorr file available for instrument\n")

        radius: float = float(config.aperture_phot_radius)
        ee_frac: float = float(config.encircled_energy_fraction)
        sky_in: float = float(config.sky_annulus_inner_radius)
        sky_out: float = float(config.sky_annulus_outer_radius)

        if ee_frac >= 0:
            radius = APPhotRoutine.radius_from_enc_energy(
                filter_string, ee_frac, ap_corr_f_name)
            if radius > 0:
                log(
                    "-> calculating aperture radius from encircled energy\n")

        if radius <= 0:
            if (radius := full_width_half_max) > 0:
                log("-> using FWHM as aperture radius\n")
            else:
                log(
                    "No valid aperture radius was detected. defaulting "
                    "to default value of 2")
                radius = 2.0

        ap_corr: float = APPhotRoutine.calc_ap_corr(
            filter_string, radius, table_f_name=ap_corr_f_name,
            verbose=config.verbose_logs)
        return radius, sky_in, sky_out, ap_corr

    # noinspection SpellCheckingInspection
    @staticmethod
    def aperture_photometry(
            detections: Table | None, image: HDUList | None,
            filter_string: str | None, log: Callable[[str], None],
            config: StarBugMainConfig, info: dict[str, str],
            full_width_half_max: float, out_dir: str | None,
            b_name: str | None, header: Header,
            background: ImageHDU | PrimaryHDU | None,
            main_image: ImageHDU | PrimaryHDU) -> Tuple[Table, ExitStates]:
        """
         Executes aperture photometry processing steps.

         Performs background subtraction and measures fluxes within specified
         aperture radii for all detected sources.

         :param detections: Table of source positions or None if empty.
         :type detections: Table | None
         :param image: Full FITS HDU list containing the data extensions.
         :type image: HDUList | None
         :param filter_string: Name of the photometric filter band used.
         :type filter_string: str | None
         :param log: Callable logging function for status tracking.
         :type log: Callable[[str], None]
         :param config: Configuration instance containing engine parameters.
         :type config: StarBugMainConfig
         :param info: Metadata dictionary detailing execution properties.
         :type info: dict[str, str]
         :param full_width_half_max: FWHM value for stellar profiles.
         :type full_width_half_max: float
         :param out_dir: Output target directory for exported tables.
         :type out_dir: str | None
         :param b_name: Base filename prefix used for exporting results.
         :type b_name: str
         :param header: FITS header containing primary data metadata.
         :type header: Header
         :param background: Estimated background map array image if loaded.
         :type background: ImageHDU | PrimaryHDU | None
         :param main_image: Primary target image array used for measurements.
         :type main_image: ImageHDU | PrimaryHDU
         :return: Extracted source table and the resulting status state.
         :rtype: Tuple[Table, ExitStates]
         """
        status: ExitStates = AperturePhotometry.check_states(
            detections, image, filter_string)
        if status is not ExitStates.EXIT_SUCCESS:
            return Table(), ExitStates.EXIT_FAIL

        assert detections is not None
        assert filter_string is not None
        assert image is not None
        detections = AperturePhotometry.remove_detection_columns(
            detections, filter_string)

        # APERTURE PHOTOMETRY
        log("\nRunning Aperture Photometry\n")

        image_data: np.ndarray
        error: np.ndarray
        mask: np.ndarray
        image_data, error, _, mask = StarBugInterface.prepare_image_arrays(
            image, log, background, header, main_image)

        # Aperture Correction
        radius, sky_in, sky_out, ap_corr = (
            AperturePhotometry.aperture_correction(
                config, info, log, filter_string, full_width_half_max))

        # Run Photometry
        app_hot: APPhotRoutine = APPhotRoutine(
            radius, sky_in, sky_out, verbose=config.verbose_logs)

        dq_flags: np.ndarray | None
        if DQ in ext_names(image):
            dq_flags = image[DQ].data.copy()
        else:
            dq_flags = None
        ap_cat: Table = app_hot(
            image_data, detections, error=error, dq_flags=dq_flags,
            ap_corr=ap_corr, sig_sky=config.sigma_sky)

        filter_string: str = filter_string if filter_string else "mag"

        # Extract magnitudes
        mag: float
        mag_err: float
        mag, mag_err = flux_to_pogson_mag(
            ap_cat[TableColumn.FLUX], ap_cat[TableColumn.E_FLUX])

        # Add columns to the catalogue
        ap_cat.add_column(Column(
            mag + config.zero_point_magnitude, filter_string))
        ap_cat.add_column(Column(
            mag_err, "e%s" % filter_string))

        # Update detections
        detections = hstack((detections, ap_cat))

        # Check for insanity
        if detections is None:
            return Table(), ExitStates.EXIT_FAIL

        if config.clean_sources:
            detections_length = len(detections)
            if (smooth_lo := config.smooth_low) != "":
                detections.remove_rows(
                    detections[TableColumn.SMOOTHNESS] < smooth_lo)
            if (smooth_hi := config.smooth_high) != "":
                detections.remove_rows(
                    detections[TableColumn.SMOOTHNESS] > smooth_hi)
            if len(detections) != detections_length:
                log("-> removing %d sources outside SMOOTH range\n"
                    % (detections_length - len(detections)))

        reindex(detections)
        detections.meta[HeaderTags.FILTER] = filter_string

        f_name = "%s/%s-ap.fits" % (out_dir, b_name)
        log("--> %s\n" % f_name)
        export_table(detections, f_name, header=header)

        return detections, ExitStates.EXIT_SUCCESS
