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
from typing import Tuple, Callable

import numpy as np
from astropy.io.fits import ImageHDU, PrimaryHDU, Header
from astropy.table import Column, Table
from astropy.wcs import WCS

from starbug2.constants import ExitStates, TableColumn
from starbug2.core.star_bug_config import StarBugMainConfig
from starbug2.routines.detection_routines import DetectionRoutine
from starbug2.utilities.filters import STAR_BUG_FILTERS
from starbug2.utilities.utils import printf, p_error


class Detect:
    """
    this class contains the main logic for running detection from
    starbug main.
    """

    @staticmethod
    def detect(
            log: Callable[[str], None],
            main_image: Callable[[], ImageHDU | PrimaryHDU],
            custom_filter: str | None, config: StarBugMainConfig,
            wcs: WCS | None,
            header: Callable[[], Header]) -> Tuple[ExitStates, Table | None]:
        """
        Full source detection routine. Saves the result as a table.

        This method executes the core star detection algorithms using
        configuration parameters provided by config, utilizing passed
        callables to dynamically fetch the target image data, logging
        mechanisms, and FITS header metadata.

        :param log: A callable function used for pipeline verbose logging.
        :type log: Callable[[str], None]
        :param main_image: A callable that retrieves the target image data
                           unit.
        :type main_image: Callable[[], ImageHDU | PrimaryHDU]
        :param custom_filter: The photometric filter band string being
                              processed.
        :type custom_filter: str | None
        :param config: Main pipeline configurations and detection thresholds.
        :type config: StarBugMainConfig
        :param wcs: World Coordinate System object for coordinate
                    translation, if available.
        :type wcs: WCS | None
        :param header: A callable that retrieves the underlying FITS header
                       metadata.
        :type header: Callable[[], Header]
        :return: A tuple containing the execution exit status state and the
                 generated source detection catalogue table (or None if
                 detection failed).
        :rtype: Tuple[ExitStates, Table | None]
        """
        log("Detecting Sources\n")
        status: ExitStates = ExitStates.EXIT_SUCCESS
        detections: Table | None = None
        assert custom_filter is not None
        if main_image():
            # noinspection SpellCheckingInspection
            detector: DetectionRoutine = DetectionRoutine(
                sig_src=config.sigma_source,
                sig_sky=config.sigma_sky,
                full_width_half_max=(
                    config.full_width_half_max_with_filter(
                        STAR_BUG_FILTERS.get(custom_filter))),
                sharp_lo=config.sharp_cutoff_low,
                sharp_hi=config.sharp_cutoff_high,
                round_1_hi=config.round1_cutoff_high,
                round_2_hi=config.round2_cutoff_high,
                smooth_lo=config.smooth_low,
                smooth_hi=config.smooth_high,
                ricker_r=config.ricker_wavelet_radius,
                do_bgd_2d=config.do_bgd_2d,
                do_con_vl=config.do_convolution,
                box_size=config.background_box_size,
                clean_src=config.clean_sources,
                verbose=config.verbose_logs)

            detections = detector(main_image().data.copy())
            assert detections is not None

            # check we have columns we need
            if not (TableColumn.X_CENTROID in detections.colnames and
                    TableColumn.Y_CENTROID in detections.colnames and
                    TableColumn.SHARPNESS in detections.colnames and
                    TableColumn.ROUNDNESS1 in detections.colnames and
                    TableColumn.ROUNDNESS2 in detections.colnames):
                printf(
                    f"Dont have the pre-requisite columns. Please ensure "
                    f"that the detections table has columns named the "
                    f"following: {TableColumn.X_CENTROID}, "
                    f"{TableColumn.Y_CENTROID}, {TableColumn.SHARPNESS}, "
                    f"{TableColumn.ROUNDNESS1}, {TableColumn.ROUNDNESS2}.")
                return ExitStates.EXIT_FAIL, detections

            # filter to just the fields we need
            detections = detections[
                TableColumn.X_CENTROID, TableColumn.Y_CENTROID,
                TableColumn.SHARPNESS, TableColumn.ROUNDNESS1,
                TableColumn.ROUNDNESS2]

            # Check for insane states
            if detections is None or wcs is None:
                return ExitStates.EXIT_FAIL, detections

            ra: np.ndarray
            dec: np.ndarray
            ra, dec = wcs.all_pix2world(
                detections[TableColumn.X_CENTROID],
                detections[TableColumn.Y_CENTROID], 0)
            detections.add_column(
                Column(ra, name=TableColumn.RA), index=2)
            detections.add_column(
                Column(dec, name=TableColumn.DEC), index=3)
            detections.meta = dict(header().items())

            # noinspection SpellCheckingInspection
            detections.meta.update({"ROUNTINE": "DETECT"})
        else:
            p_error("Something went wrong.\n")
            status = ExitStates.EXIT_FAIL
        return status, detections
