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
from typing import List, Optional, Union

import numpy as np
from matplotlib.axis import Axis
from photutils.background import Background2D, BackgroundBase
from photutils.aperture import CircularAperture, ApertureMask
from astropy.table import Table
from starbug2.constants import TableColumn
from starbug2.utilities.utils import Loading, printf, warn


class BackGroundEstimateRoutine(BackgroundBase):
    def __init__(
            self, source_list: Table | None,
            box_size: int = 2, full_width_half_max: float = 2.0,
            sig_sky: float = 2.0, bgd_r: float = -1.0,
            profile_scale: float = 1.0, profile_slope: float = 0.5,
            verbose: bool | int = 0,
            bgd: Optional[Background2D] = None) -> None:
        """
        Diffuse background emission estimator run by starbug.

        :param source_list: List of sources in the image in a table
                            containing X_CENTROID Y_CENTROID
        :type source_list: astropy.table.Table or None
        :param box_size: Size of the kernel to pass over the image in
                         un-sharp masking
        :type box_size: int
        :param full_width_half_max: Source full width half maximum in the
                                    image.
        :type full_width_half_max: float.
        :param sig_sky: Sigma threshold for background clipping.
        :type sig_sky: float.
        :param bgd_r: Fixed aperture mask radius around each source.
        :type bgd_r: float.
        :param profile_scale: Scaling factor for the aperture mask radius
                              profile.
        :type profile_scale: float.
        :param profile_slope: Slope of the aperture mask radius profile.
        :type profile_slope: float.
        :param verbose: Set whether to print verbose output information.
        :type verbose: bool.
        """
        self._source_list: Table | None = source_list
        self._box_size: int = box_size
        self._full_width_half_max: float = full_width_half_max
        self._sig_sky: float = sig_sky
        self._bgd_r: float = bgd_r
        self._a: float = profile_scale
        self._b: float = profile_slope
        self._verbose: int | bool = verbose
        self._background: Background2D | None = bgd
        super().__init__()

    def calc_peaks(self, im: np.ndarray) -> np.ndarray:
        """
        Determine peak pixel value for each source in xy
        :param im: ??????
        :return: peaks
        :rtype: np.array
        """

        assert self._source_list is not None
        x: Table = self._source_list[TableColumn.X_CENTROID]
        y: Table = self._source_list[TableColumn.Y_CENTROID]
        apertures: List[ApertureMask] = CircularAperture(
            np.array((x, y)).T, 2).to_mask()
        peaks: np.ndarray = np.full(len(x), np.nan)

        i: int
        mask: ApertureMask
        for i, mask in enumerate(apertures):
            peaks[i] = np.nanmax(mask.multiply(im))
        return peaks

    def log(self, msg: str) -> None:
        """
        log this message.
        :param msg: the message to log
        :return: None
        """
        if self._verbose:
            printf(msg)

    def __call__(
            self, data: np.ndarray | None,
            axis: Optional[Axis] = None, masked: bool = False,
            output: Optional[str] = None) -> Background2D | None:
        """
        does background estimation routine.

        :param data: the data to process
        :type data: np.ndarray
        :param axis: the axis
        :type axis: matplotlib.axis.Axis or None
        :param masked: if the data is masked
        :type masked: bool
        :param output: if the data should be outputted
        :type output: str
        :return: the new background 2D object.
        :rtype: Background2D
        """
        if self._source_list is None or data is None:
            return self._background
        _data: np.ndarray = np.copy(data)

        x_grid: np.ndarray
        y_grid: np.ndarray
        x_grid, y_grid = np.ogrid[:data.shape[1], :data.shape[0]]

        default_r: float = 2 * self._full_width_half_max

        rlist: np.ndarray
        if self._bgd_r and self._bgd_r > 0:
            self.log(
                "-> using BGD_R=%g masking aperture radii\n" % self._bgd_r)
            rlist = self._bgd_r * np.ones(len(self._source_list))
        else:
            if TableColumn.FLUX in self._source_list.colnames:
                self.log("-> calculating source aperture mask radii\n")
                sky: Union[np.ndarray, float] = (
                    self._source_list["sky"]
                    if "sky" in self._source_list.colnames else 1.0)
                rlist = (
                    self._a * self._full_width_half_max * (
                        np.log(self._source_list[TableColumn.FLUX] / sky))
                    ** self._b)
                rlist[np.isnan(rlist)] = default_r

                if output:
                    with open(output, 'w') as fp:
                        for i in range(len(rlist)):
                            radius_val: float = float(rlist[i])
                            x_cen: float = float(
                                self._source_list[i][TableColumn.X_CENTROID])
                            y_cen: float = float(
                                self._source_list[i][TableColumn.Y_CENTROID])

                            fp.write(
                                "circle %f %f %f #color=green;" % (
                                    1 + x_cen, 1 + y_cen, radius_val))
                            fp.write(
                                "annulus %f %f %f %f #color=white;" % (
                                    1 + x_cen, 1 + y_cen, 1.5 * radius_val,
                                    1.5 * radius_val + 1))
                    self.log("-> exporting check file \"%s\"\n" % output)
            else:
                warn("Unable to calculate aperture mask sizes, "
                     "add '-A' to starbug command.\n")
                rlist = default_r * np.ones(len(self._source_list))

        dimension: int = 50
        load: Loading = Loading(
            len(self._source_list), msg="masking sources", res=10)

        r: float
        src: Table
        for r, src in zip(rlist, self._source_list):  # type: ignore

            rin: float = 1.5 * r
            rout: float = rin + 1

            x: int = int(round(src[TableColumn.X_CENTROID]))
            y: int = int(round(src[TableColumn.Y_CENTROID]))
            _X: np.ndarray = x_grid[
                max(x - dimension, 0): min(x + dimension, data.shape[1])]
            _Y: np.ndarray = y_grid[
                :, max(y - dimension, 0): min(y + dimension, data.shape[0])]

            radius: np.ndarray = np.sqrt(
                (_X - src[TableColumn.X_CENTROID]) ** 2 +
                (_Y - src[TableColumn.Y_CENTROID]) ** 2)

            mask: np.ndarray = (radius < r)
            annuli_mask: np.ndarray = ((radius > rin) & (radius < rout))

            tmp: np.ndarray = _data[_Y, _X]
            tmp[mask] = np.median(data[_Y, _X][annuli_mask])
            _data[_Y, _X] = tmp

            load()

            # This will slow the thing down quite a lot
            if self._verbose:
                load.show()
        if self._verbose:
            printf("-> estimating bgd2d\n")
        self._background = Background2D(_data, self._box_size)
        return self._background

    def calc_background(
            self, data: np.ndarray,
            axis: Optional[int] = None,
            masked: Optional[bool] = None) -> np.ndarray:
        """
        Calculate the background value.

        Parameters
        --------------
        data : array_like or `~numpy.ma.MaskedArray`
            The array for which to calculate the background value.

        axis : int or `None`, optional
            The array axis along which the background is calculated. If
            `None`, then the entire array is used.

        masked : bool, optional
            If `True`, then a `~numpy.ma.MaskedArray` is returned. If
            `False`, then a `~numpy.ndarray` is returned, where masked
            values have a value of NaN. The default is `False`.

        Returns
        -----------
        result : float, `~numpy.ndarray`, or `~numpy.ma.MaskedArray`
            The calculated background value. If ``masked`` is
            `False`, then a `~numpy.ndarray` is returned, otherwise a
            `~numpy.ma.MaskedArray` is returned. A scalar result is
            always returned as a float.
        """
        if self._background is None:
            self.__call__(data)
        if self._background is None:
            raise Exception(
                "the background file didnt get created for some reason")
        return self._background.background
