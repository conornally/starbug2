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
import sys

import numpy as np
from astropy.table import Column, hstack, Table, QTable
from photutils.aperture import CircularAperture, aperture_photometry
from photutils.psf import PSFPhotometry, SourceGrouper, ImagePSF

from starbug2.constants import TableColumn
from starbug2.utilities.utils import printf, p_error, warn


class _Grouper(SourceGrouper):
    """
    Overloaded SourceGrouper. This class gives a starbug warning into stderr
     if there are more than CRITICAL_VAL source in any given group. Then
     returns the original __call__ function results.

    Parameters:
    -----------
    min_separation : float > 0
        The minimum distance (in pixels) such that any two sources
        separated by less than this distance will be placed in the same
        group if the ``min_size`` criteria is also met.
    """
    CRITICAL_VAL = 25
    min_separation = 0

    def __init__(self, min_separation: float) -> None:
        super().__init__(min_separation)

    def __call__(self, x: np.ndarray, y: np.ndarray,
                 return_groups_object: bool = False) -> np.ndarray:
        res: np.ndarray = super().__call__(x, y)
        n: int
        if n := sum(np.bincount(res) > self.CRITICAL_VAL):  # noqa
            warn("Source grouper has %d groups larger than %d. Consider"
                 " reducing \"CRIT_SEP=%g\" or fitting might take a long"
                 " time.\n" % (n, self.CRITICAL_VAL, self.min_separation))
        return res


class PSFPhotRoutine(PSFPhotometry):
    def __init__(
            self, psf_model: ImagePSF,
            fit_shape: int | tuple[int, int],
            app_hot_r: float = 3.0,
            min_separation: float = 8.0,
            force_fit: int | bool = False,
            background: np.ndarray | None = None,
            verbose: int | bool = 1) -> None:
        # noinspection SpellCheckingInspection
        """
        PSF Photometry routine called by starbug

        :param psf_model: Model PSF to be used in the fitting
        :type psf_model: photutils.psf.FittableImageModel
        :param fit_shape: Size of PSF to use in pixels. Must be less than or
                          equal to the size of psf_model
        :type fit_shape: int or tuple
        :param min_separation: Minimum source separation for source grouper
                               (pixels)
        :type min_separation: float
        :param app_hot_r: Aperture radius to be used in initial guess
                          photometry.
        :type app_hot_r: float.
        :param force_fit: Conduct forced centroid PSF fitting.
        :type force_fit: bool or int (with values 0 or 1).
        :param background: 2D array with the same dimensions as the data used
                           in fitting.
        :type background: numpy.ndarray.
        :param verbose: Show verbose outputs.
        :type verbose: bool or int
        """
        self._verbose: int | bool = verbose
        self._force_fit: bool | int = force_fit
        self._background: np.ndarray | None = background

        grouper: _Grouper = _Grouper(min_separation)

        if force_fit:
            psf_model.x_0.fixed = True
            psf_model.y_0.fixed = True

        super().__init__(
            psf_model=psf_model, fit_shape=fit_shape, finder=None,
            progress_bar=verbose, aperture_radius=app_hot_r, grouper=grouper)

        if self._verbose:
            printf("-> source group separation: %g\n" % min_separation)

    def __call__(
            self, image: np.ndarray,
            init_params: Table | None = None,
            error: np.ndarray | None = None,
            mask: np.ndarray | None = None):
        """
        runs the psf phot routine.
        :param image: the image to process.
        :type image: np.ndarray
        :param init_params: the init params.
        :type init_params: Table
        :param error: the error.
        :type error: np.array
        :param mask: the mask.
        :type mask: np.array
        :return: the processed table.
        :rtype:  astropy.table.Table
        """
        return self.do_photometry(image, init_params, error, mask)

    def do_photometry(
            self,
            image: np.ndarray,
            init_params: Table | None = None,
            error: np.ndarray | None = None,
            mask: np.ndarray | None = None) -> Table | None:
        """
        does the photometry
        :param image: the image to process.
        :type image: np.ndarray
        :param init_params: the init params.
        :type init_params: Table
        :param error: the error.
        :type error: np.array
        :param mask: the mask.
        :type mask: np.array
        :return: the processed table.
        :rtype:  astropy.table.Table
        """

        if (init_params is None or len(init_params) == 0
                or mask is None or error is None):
            p_error("Must include source list and a mask\n")
            return None

        # Removing completely masked sources
        apertures: CircularAperture = CircularAperture(
            [(row[TableColumn.X_INIT],
              row[TableColumn.Y_INIT]) for row in init_params],
            self.aperture_radius)
        ap_masks: QTable = aperture_photometry(~mask, apertures)
        init_params.remove_rows(ap_masks["aperture_sum"] == 0)

        # bad errors should be big not small
        error[error == 0] = sys.maxsize

        if self._background is not None:
            image = image - self._background
        if self._verbose:
            printf("-> fitting %d sources\n" % len(init_params))
        cat: QTable = super().__call__(
            image, mask=mask, init_params=init_params, error=error)

        d: np.ndarray = np.sqrt((
            (cat[TableColumn.X_INIT] - cat[TableColumn.X_FIT]) ** 2.0 +
            (cat[TableColumn.Y_INIT] - cat[TableColumn.Y_FIT]) ** 2.0))

        # noinspection SpellCheckingInspection
        cat.add_column(Column(d, name=TableColumn.XY_DEV))

        if TableColumn.FLUX_ERR not in cat.colnames:
            cat.add_column(
                Column(np.full(len(cat), np.nan), name=TableColumn.E_FLUX))
            warn("Something went wrong with PSF error fitting\n")
        else:
            cat.rename_column(TableColumn.FLUX_ERR, TableColumn.E_FLUX)

        cat.rename_column(TableColumn.FLUX_FIT, TableColumn.FLUX)

        # noinspection SpellCheckingInspection
        keep: list[str] = [
            TableColumn.X_FIT, TableColumn.Y_FIT, TableColumn.FLUX,
            TableColumn.E_FLUX, TableColumn.XY_DEV, TableColumn.Q_FIT]
        return hstack((init_params, cat[keep]))
