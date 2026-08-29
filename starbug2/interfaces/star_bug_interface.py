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

from abc import ABC, abstractmethod
from typing import Tuple, Callable

from astropy.table import Table
from astropy.io.fits import PrimaryHDU, ImageHDU, HDUList, Header

import numpy as np

from starbug2.constants import (
    ExitStates, ImageHeaderTags, AREA, ERR, DQ, DQFlags)
from starbug2.utilities.utils import get_mj_ysr2jy_scale_factor, ext_names


class StarBugInterface(ABC):

    @abstractmethod
    def log(self, msg: str) -> None:
        """
        Print message if in verbose mode

        :param msg: Message to print out
        :type msg: str
        :return: None
        """
        pass

    @abstractmethod
    def load_image(self, f_name: str) -> None:
        """
        Given f_name, load the image into starbug to be worked on.

        :param f_name: Filename of fits image (with any number of extensions).
            If using a non-standard HDU index, set the name or index of the
            extension with "HDU_NAME=XXX" in the parameter file.
        :type f_name: str
        :return: None
        """
        pass

    @abstractmethod
    def load_ap_file(self, f_name: str | None = None) -> None:
        """
        Load an AP_FILE to be used during photometry

        :param f_name: Filename for fits table containing source coordinates.
         These coordinates can be x-centroid / y-centroid, x_init / y_init,
          x_0, y_0 or RA/DEC. The latter is used if starbug gets "USE_WCS=1"
          in the parameter file.
        :type f_name: str or None
        :return: None
        """
        pass

    @abstractmethod
    def load_bgd_file(self, f_name: str | None = None) -> None:
        """
        Load a BGD_FILE to be used during photometry

        :param f_name: Filename of fits image the same dimensions as the
                       main image
        :type f_name: str or None
        :return: None
        """
        pass

    @abstractmethod
    def detect(self) -> ExitStates:
        """
        Full source detection routine. Saves the result as a table
        self._detections

        :return: The execution status (0 for success, non-zero for failure)
        :rtype: ExitStates
        """
        pass

    @abstractmethod
    def aperture_photometry(self) -> ExitStates:
        """
        Executes aperture photometry

        :return: 0 for success, 1 for failure
        :rtype: ExitStates
        """
        pass

    @abstractmethod
    def bgd_estimate(self) -> ExitStates:
        """
        Estimate the background of the active image
        Saves the result as an ImageHDU self._background

        :return: The execution status (0 for success, non-zero for failure)
        :rtype: ExitStates
        """
        pass

    @abstractmethod
    def bgd_subtraction(self) -> ExitStates:
        """
        Internally subtract a background array from an image array

        :return: 0 for success, 1 otherwise
        :rtype: ExitStates
        """
        pass

    @abstractmethod
    def psf_photometry_routine(self) -> ExitStates:
        """
        Full photometry routine
        Saves the result as a table self._psf_catalogue,
        Additionally it appends a residual Image onto the
        self._residuals HDUList

        :return: 0 for success, 1 otherwise
        :rtype: ExitStates
        """
        pass

    @abstractmethod
    def source_geometry(self) -> None:
        """
        Calculate source geometry stats for a given image and source list

        :return: None
        """
        pass

    @abstractmethod
    def verify(self) -> ExitStates:
        """
        This simple function verifies that everything necessary has been
        loaded properly

        :return: 0 on success, 1 on failure
        :rtype: ExitStates
        """
        pass

    @abstractmethod
    def header(self) -> Header:
        """
        Construct relevant base header information for routine products

        :return: Header file containing a series of relevant information
        :rtype: Header
        """
        pass

    @property
    @abstractmethod
    def info(self) -> dict[str, str]:
        """
        Get some useful information from the image header file.

        :return: Extracted keys and elements from the image header.
        :rtype: dict of str to str
        """
        pass

    @abstractmethod
    def main_image(self) -> ImageHDU | PrimaryHDU:
        # noinspection SpellCheckingInspection
        """
        Automagically find the main image array to use
        Order of importance is:
        > self._nHDU (if set)
        > param[ HDUNAME ]
        > SCI, BGD, RES
        > first ImageHDU
        > image[0]

        :return: The main image array extension
        :rtype: ImageHDU or PrimaryHDU
        """
        pass

    @property
    @abstractmethod
    def filter(self) -> str | None:
        pass

    @property
    @abstractmethod
    def n_hdu(self) -> int:
        pass

    @property
    @abstractmethod
    def image(self) -> HDUList | None:
        pass

    @image.setter
    @abstractmethod
    def image(self, new_image: HDUList) -> None:
        pass

    @property
    @abstractmethod
    def psf_catalogue(self) -> Table | None:
        pass

    @property
    @abstractmethod
    def psf(self) -> np.ndarray | None:
        pass

    @psf.setter
    @abstractmethod
    def psf(self, new_value: np.ndarray) -> None:
        pass

    @property
    @abstractmethod
    def f_name(self) -> str | None:
        pass

    @property
    @abstractmethod
    def detections(self) -> Table | None:
        pass

    @detections.setter
    @abstractmethod
    def detections(self, new_detections: Table) -> None:
        pass

    @property
    @abstractmethod
    def out_dir(self) -> str | None:
        pass

    @staticmethod
    def prepare_image_arrays(
            image: HDUList | None, log: Callable[[str], None],
            background: ImageHDU | PrimaryHDU | None,
            header: Header, main_image: ImageHDU | PrimaryHDU) -> (
            Tuple[np.ndarray, np.ndarray, np.ndarray | None, np.ndarray]):
        """
        Make a copy of the original image, and prepare the other image arrays

        Extracts and copies primary data, handles associated variance or error
        extensions, validates background arrays, and generates pixel masks.

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
        :return: A tuple containing (image, error, bgd, mask) arrays.
        :rtype: Tuple[np.ndarray, np.ndarray, np.ndarray | None, np.ndarray]
        """
        # Collect scale factor
        scale_factor: int | float
        if header.get(ImageHeaderTags.BUN_IT) == "MJy/sr":
            scale_factor = get_mj_ysr2jy_scale_factor(main_image)
            log(
                "-> converting unit from MJy/sr to Jr with factor: %e\n"
                % scale_factor)
        else:
            scale_factor = 1

        image_data: np.ndarray = main_image.data.copy() * scale_factor

        # Scale by area
        extension_names: list[str] = ext_names(image)
        assert image is not None
        if AREA in extension_names:
            # AREA distortion correction
            image_data *= image[AREA].data

        # Collect and scale error
        error: np.ndarray
        if ERR in extension_names and np.shape(image[ERR]):
            error = image[ERR].data.copy() * scale_factor
        else:
            error = np.sqrt(np.abs(image_data))

        # Create mask
        mask: np.ndarray
        if DQ in extension_names:
            mask = (
                image[DQ].data
                & (DQFlags.DQ_DO_NOT_USE | DQFlags.DQ_SATURATED))
            mask = mask.astype(bool)
        else:
            mask = (np.isnan(image_data) | np.isnan(error))

        # Collect and scale background array
        bgd: np.ndarray | None
        if background is not None:
            bgd = background.data.copy() * scale_factor
        else:
            bgd = None

        return image_data, error, bgd, mask
