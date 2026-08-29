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

from typing import Optional, cast
from collections.abc import Callable

import numpy as np
from scipy.ndimage import convolve
from skimage.feature import match_template

from astropy.stats import sigma_clipped_stats
from astropy.coordinates import SkyCoord
from astropy.table import Column, Table, vstack
from astropy.convolution import RickerWavelet2DKernel
from astropy.units import Quantity

from photutils.background import Background2D
from photutils.detection import StarFinderBase, DAOStarFinder, find_peaks

from starbug2.constants import TableColumn
from starbug2.routines.source_properties import SourceProperties
from starbug2.utilities.utils import printf


class DetectionRoutine(StarFinderBase):
    def __init__(
            self, sig_src: float = 5.0, sig_sky: float = 3.0,
            full_width_half_max: float = 2.0, sharp_lo: float = 0.2,
            sharp_hi: float = 1, round_1_hi: float = 1, round_2_hi: float = 1,
            smooth_lo: float = -np.inf, smooth_hi: float = np.inf,
            ricker_r: float = 1.0, verbose: int | bool = 0,
            clean_src: bool | int = 1, do_bgd_2d: bool | int = 1,
            box_size: int = 2, do_con_vl: bool | int = 1) -> None:
        # noinspection SpellCheckingInspection
        """
        Detection routine

        A standalone detection that runs on a 2D image.
        It uses DAOStarFinder as the base for peak detection but run
        several times on a series of background subtracted images.
        Each run the background subtraction is different, bringing out a
        different set of sources

        :param sig_src:  The detection flux threshold, this sets the number of
                         sigma above the  image median flux that a source must
                         be brighter than. Default:  sig_src=5 is a "solid"
                         detection.
        :type sig_src: float
        :param sig_sky: The number of sigma above the image median flux that
                        is still considered "sky". Pixels below this will be
                        cut out during the  detection steps.
        :type sig_sky: float
        :param full_width_half_max: Full width half maximum of a standard
                                    source in the image.
        :type full_width_half_max: float
        :param sharp_lo: Lowest bound for a source "sharpness".
        :type sharp_lo: float
        :param sharp_hi: Upper bound for a source "sharpness".
        :type sharp_hi: float
        :param round_1_hi: Upper bound for a source "roundness1", this
                           distribution is symmetric and the lower bound is
                           taken as negative round_1_hi.
        :type round_1_hi: float
        :param round_2_hi: Upper bound for a source "roundness2", this
                           distribution is symmetric and the lower bound is
                           taken as negative round_2_hi.
        :type round_2_hi: float
        :param smooth_lo: Lower bound for source "smoothness".
        :type smooth_lo: float
        :param smooth_hi: Upper bound for source "smoothness".
        :type smooth_hi: float
        :param ricker_r: Pixel radius for the wavelet used in the CONVL
                         detection step.
        :type ricker_r: float
        :param verbose: Set whether to print verbose output information.
        :type verbose: bool or int
        :param clean_src: Set whether to "clean" the catalogue after detection
                          based on the above source geometric properties.
        :type clean_src: bool
        :param do_bgd_2d: Set whether to run the BGD2D detection step.
        :type do_bgd_2d: bool
        :param box_size: Set kernel size for BGD2D background measuring step.
        :type box_size: int
        :param do_con_vl: Set whether to run the CONVL detection step.
        :type do_con_vl: bool
        """
        self.sig_src: float = sig_src
        self.sig_sky: float = sig_sky
        self.full_width_half_max: float = full_width_half_max
        self.sharp_hi: float = sharp_hi
        self.sharp_lo: float = sharp_lo
        self.round_1_hi: float = (
            round_1_hi if round_1_hi is not None else np.inf)
        self.round_2_hi: float = (
            round_2_hi if round_2_hi is not None else np.inf)
        self.smooth_lo: float = (
            smooth_lo if smooth_lo is not None else -np.inf)
        self.smooth_hi: float = (
            smooth_hi if smooth_hi is not None else np.inf)

        self.ricker_r: float = ricker_r
        self.clean_src: bool | int = clean_src

        self.catalogue: Table = Table()
        self.verbose: bool | int = verbose

        self.do_bgd_2d: bool | int = do_bgd_2d
        self.box_size: int = box_size
        self.do_con_vl: bool | int = do_con_vl

    def detect(self, data: np.ndarray,
               bkg_estimator: Optional[Callable[
                   [np.ndarray], np.ndarray]] = None,
               xy_coords: Table | None = None,
               use_find_peaks: bool = False) -> Table:
        """
        The core detection step (DAOStarFinder)

        :param data: Image array to detect on
        :type data: numpy.ndarray or array.pyi
        :param bkg_estimator: Function to call to generate the background
                              array the same shape as data array
        :type bkg_estimator:  callable function
        :param xy_coords: Table of initial guesses (x_centroid, y_centroid)
        :type xy_coords: `astropy.table.Table`
        :param use_find_peaks: if true, use the photutils findpeaks method
                               else use the DAOStarFinder method
        :type use_find_peaks: bool
        :return: Source list Table
        :rtype: astropy.Table
        """
        bkg: np.ndarray = np.zeros(data.shape)
        if bkg_estimator:
            bkg = bkg_estimator(data)

        median_stat: float
        std: float
        catalogue: Table | None
        _, median_stat, std = sigma_clipped_stats(data, sigma=self.sig_sky)
        if use_find_peaks:
            catalogue = find_peaks(
                data - bkg, median_stat + std * self.sig_src, box_size=11)

        else:
            round_hi: float = max((self.round_1_hi, self.round_2_hi))
            find: DAOStarFinder = DAOStarFinder(
                threshold=std * self.sig_src, fwhm=self.full_width_half_max,
                sharplo=self.sharp_lo, sharphi=self.sharp_hi,
                roundlo=-round_hi, roundhi=round_hi, peak_max=np.inf,
                xycoords=xy_coords)
            catalogue = find.find_stars(data - bkg)

        # if no results, create a table reflecting no results. with expected
        # column names.
        if catalogue is None or len(catalogue) == 0:
            col_names = [
                TableColumn.ID, TableColumn.X_PEAK, TableColumn.Y_PEAK,
                TableColumn.X_CENTROID, TableColumn.Y_CENTROID]
            col_data = [list(), list(), list(), list(), list()]
            catalogue: Table = Table(col_data, names=col_names)
        return catalogue

    def bkg2d(self, data: np.ndarray) -> np.ndarray:
        """
        Calculates a 2D background array map.

        :param data: the data to apply background 2d to.
        :return: background
        :rtype: numpy.ndarray
        """
        return Background2D(data, self.box_size, filter_size=3).background

    def match(self, base: Table, cat: Table) -> Table:
        """
         Internal function to class
        Used to match detections from separate background subtracted images
        into the main catalogue. This will append a source if its matched
        separation  is above the threshold = self.full_width_half_max

        :param base: Base catalogue to match to.
        :type base: astropy.table.Table
        :param cat: Catalogue to be matched
        :type cat: astropy.table.Table
        :return: The matched catalogue
        :rtype: astropy.Table
        """
        base_sky: SkyCoord = SkyCoord(
            x=base[TableColumn.X_CENTROID],
            y=base[TableColumn.Y_CENTROID],
            z=np.zeros(len(base)),
            representation_type="cartesian")
        cat_sky: SkyCoord = SkyCoord(
            x=cat[TableColumn.X_CENTROID],
            y=cat[TableColumn.Y_CENTROID],
            z=np.zeros(len(cat)),
            representation_type="cartesian")

        dist: Quantity
        _, _, dist = cat_sky.match_to_catalog_3d(base_sky)
        mask: np.ndarray = dist.to_value() > self.full_width_half_max
        return vstack((base, cat[mask]))

    def _clean_up_data(
            self, data: np.ndarray | None,
            mask: Optional[np.ndarray] = None) -> np.ndarray | None:
        """
        sets up data for detections
        :param data: 2D image array to detect on
        :type data: numpy.ndarray or None
        :param mask: Pixels to mask out on the data array
        :type mask: numpy.ndarray
        :return: the masked data
        """
        if data is None:
            return None
        if mask is None:
            mask = np.where(np.isnan(data))

        median_stat: float
        _, median_stat, _ = sigma_clipped_stats(data, sigma=self.sig_sky)
        data[mask] = median_stat
        return data

    def _do_plain_detections(self, data: np.ndarray) -> None:
        """
        do plain detections.

        :param data: 2D image array to detect on
        :type data: numpy.ndarray or None
        :return: None
        """
        self.catalogue = self.detect(data)
        if self.verbose:
            printf("-> [PLAIN] pass: %d sources\n" % len(self.catalogue))

    def _do_background_detection(self, data: np.ndarray) -> None:
        """
        executes the background detection.

        :param data: 2D image array to detect on
        :type data: numpy.ndarray or None
        :return: None
        """

        # execute background search and then feed it to a match with
        # plain results.
        background_results: Table = self.detect(data, self.bkg2d)
        if len(self.catalogue) == 0:
            self.catalogue = background_results
        else:
            self.catalogue = self.match(self.catalogue, background_results)
        if self.verbose:
            printf("-> [BGD2D] pass: %d sources\n" % len(self.catalogue))

    def _do_differential_detection(self, data: np.ndarray) -> None:
        """
        executes the 2nd differential detection.

        :param data: 2D image array to detect on
        :type data: numpy.ndarray or None
        :return: None
        """

        kernel: RickerWavelet2DKernel = (
            RickerWavelet2DKernel(self.ricker_r))
        conv: np.ndarray = cast(np.ndarray, convolve(data, kernel.array))
        corr: np.ndarray = cast(np.ndarray, match_template(
            conv / np.amax(conv), kernel.array))
        detections: Table = self.detect(corr, use_find_peaks=True)
        if detections:
            detections[TableColumn.X_PEAK] += kernel.shape[0] // 2
            detections[TableColumn.Y_PEAK] += kernel.shape[0] // 2
            detections.rename_columns(
                (TableColumn.X_PEAK, TableColumn.Y_PEAK),
                (TableColumn.X_CENTROID, TableColumn.Y_CENTROID))
            if len(self.catalogue) == 0:
                self.catalogue = detections
            else:
                self.catalogue = self.match(self.catalogue, detections)
        if self.verbose:
            # noinspection SpellCheckingInspection
            printf("-> [CONVL] pass: %d sources\n" % len(self.catalogue))

    def _update_mask_for_cleaning_src(self, mask: np.ndarray) -> np.ndarray:
        if (TableColumn.SHARPNESS in self.catalogue.colnames and
                TableColumn.ROUNDNESS1 in self.catalogue.colnames and
                TableColumn.ROUNDNESS2 in self.catalogue.colnames):
            mask &= (
                (self.catalogue[TableColumn.SHARPNESS] > self.sharp_lo)
                & (self.catalogue[TableColumn.SHARPNESS] < self.sharp_hi)
                & (self.catalogue[TableColumn.ROUNDNESS1] > -self.round_1_hi)
                & (self.catalogue[TableColumn.ROUNDNESS1] < self.round_1_hi)
                & (self.catalogue[TableColumn.ROUNDNESS2] > -self.round_2_hi)
                & (self.catalogue[TableColumn.ROUNDNESS2] < self.round_2_hi))
        else:
            printf(f"Incorrect columns for cleaning source. not doing so. "
                   f"Please ensure the catalogue table contains the following"
                   f" column names: [{TableColumn.SHARPNESS}, "
                   f"{TableColumn.ROUNDNESS1}, {TableColumn.ROUNDNESS2}]")
        return mask

    def _reshape_detections(self, data: np.ndarray) -> None:
        """
        reshapes the sharp and round at detected locations.

        :param data: 2D image array to detect on
        :type data: numpy.ndarray or None
        :return: None
        """
        tmp: Table | None = (
            SourceProperties(data, self.catalogue, verbose=self.verbose)
            .calculate_geometry(self.full_width_half_max))
        if tmp is not None:
            self.catalogue = tmp

        mask: np.ndarray = (
            ~np.isnan(self.catalogue[TableColumn.X_CENTROID]) &
            ~np.isnan(self.catalogue[TableColumn.Y_CENTROID]))

        if self.clean_src:
            self._update_mask_for_cleaning_src(mask)

        if self.verbose:
            printf("-> cleaning %d unlikely point sources\n" % sum(~mask))
        self.catalogue = self.catalogue[mask]

        if self.verbose:
            printf("Total: %d sources\n" % len(self.catalogue))

        self.catalogue.replace_column(
            "id", Column(range(1, 1 + len(self.catalogue))))

    def find_stars(
            self, data: np.ndarray | None,
            mask: Optional[np.ndarray] = None) -> Table:
        """
        This routine runs source detection several times, but on a different
        form of the data array each time. Each form has been "skewed" somehow
        to brighten the most faint sources and flatten the differential
        background.

        1:Plain detections
        2:Subtract Background estimation
        3:RickerWave convolution

        :param data: 2D image array to detect on
        :type data: numpy.ndarray or None
        :param mask: Pixels to mask out on the data array
        :type mask: numpy.ndarray
        :return: the catalogue containing stars.
        :rtype: astropy.Table
        """
        data: np.ndarray | None = self._clean_up_data(data, mask)
        if data is None:
            printf("No data, not trying to detect")
            return Table()

        self._do_plain_detections(data)
        if self.do_bgd_2d:
            self._do_background_detection(data)

        # 2nd order differential detection
        if self.do_con_vl:
            self._do_differential_detection(data)

        # Now with xy-coords DAOStarfinder will refit the sharp and round
        # values at the detected locations
        self._reshape_detections(data)
        return self.catalogue
