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
from typing import Tuple

import numpy as np

from astropy.stats import sigma_clip
from astropy.table import Column, Table, Row, QTable

from photutils.aperture import (
    CircularAperture, CircularAnnulus, aperture_photometry, ApertureMask)

from starbug2.constants import (
    SourceFlags, DQFlags, TableColumn, HeaderTags, Modes, QTableColNames)
from starbug2.utilities.utils import printf, p_error, warn


class APPhotRoutine:

    @staticmethod
    def calc_ap_corr(
            filter_string: str, radius: float,
            table_f_name: str | None = None,
            verbose: int | bool = 0) -> float:
        """
        Using CRDS ap_corr table, fit a curve to the radius vs ap_corr
        columns and then return ap_corr to respective input radius

        :param filter_string: the filter string
        :type filter_string: str
        :param radius: the radius
        :type radius: float
        :param table_f_name: the table file name
        :type table_f_name: str
        :param verbose: int for verbose.
        :type verbose: int
        :return: the ap_corr float.
        :rtype: float
        :raises FileNotFoundError when the table_f_name does not exist
        """

        if not table_f_name or not os.path.exists(table_f_name):
            raise FileNotFoundError("cant find the table filename")

        tmp: Table = Table.read(table_f_name, format="fits")
        t_ap_corr: Table
        if HeaderTags.FILTER_LOWER in tmp.colnames:
            t_ap_corr = tmp[(tmp[HeaderTags.FILTER_LOWER] == filter_string)]
        else:
            t_ap_corr = tmp

        if TableColumn.PUPIL in t_ap_corr.colnames:
            t_ap_corr = t_ap_corr[t_ap_corr[TableColumn.PUPIL] == Modes.CLEAR]

        ap_corr: float = float(np.interp(
            radius, t_ap_corr[TableColumn.RADIUS],
            t_ap_corr[TableColumn.AP_CORR]))
        if verbose:
            printf("-> estimating aperture correction: %.3g\n" % ap_corr)
        return ap_corr

    @staticmethod
    def ap_corr_from_enc_energy(
            filter_string: str, encircled_energy: np.ndarray,
            table_f_name: str | None = None, verbose: int | bool = 0) -> (
                Tuple[np.ndarray, np.ndarray]):
        """
        Rather than fitting radius to the AP_CORR CRDS, use the closes
        Encircled energy value

        :param filter_string: the filter string
        :type filter_string: str
        :param encircled_energy: the encircled energy
        :type encircled_energy: np.array
        :param table_f_name: the table file name
        :type table_f_name: str
        :param verbose: int for verbose.
        :type verbose: int | bool
        :return: tuple of ap_corr and radius or ExitFail.
        :rtype: tuple of np.array and np.array or int.
        :raises: FileNotFoundError when the table_f_name does not exist.
        """
        if not table_f_name or not os.path.exists(table_f_name):
            raise FileNotFoundError("cannot find table f name")

        tmp: Table = Table.read(table_f_name, format="fits")

        if HeaderTags.FILTER_LOWER in tmp.colnames:
            t_ap_corr = tmp[(tmp[HeaderTags.FILTER_LOWER] == filter_string)]
        else:
            t_ap_corr = tmp

        line: Row = t_ap_corr[(np.abs(
            t_ap_corr[TableColumn.EE_FRACTION] - encircled_energy)).argmin()]
        if verbose:
            printf(
                "-> best matching encircled energy %.1f,"
                " with radius %g pixels\n" % (
                    line[TableColumn.EE_FRACTION], line[TableColumn.RADIUS]))
            printf(
                "-> using aperture correction: %f\n" %
                line[TableColumn.AP_CORR])

        return line[TableColumn.AP_CORR], line[TableColumn.RADIUS]

    @staticmethod
    def radius_from_enc_energy(
            filter_string: str, ee_frac: float,
            table_f_name: str | None) -> float:
        """
        """
        if not table_f_name or not os.path.exists(table_f_name):
            raise FileNotFoundError("cannot find table f name")
        t_ap_corr = Table.read(table_f_name, format="fits")

        if (len({TableColumn.EE_FRACTION, TableColumn.RADIUS} &
                set(t_ap_corr.col_names)) != 2):
            raise Exception("invalid col_names size.")

        # Crop down table
        if HeaderTags.FILTER_LOWER in t_ap_corr.col_names:
            t_ap_corr = t_ap_corr[
                (t_ap_corr[HeaderTags.FILTER_LOWER] == filter_string)]

        # Crop down table
        if TableColumn.PUPIL in t_ap_corr.col_names:
            t_ap_corr = t_ap_corr[t_ap_corr[TableColumn.PUPIL] == Modes.CLEAR]

        return float(
            np.interp(
                ee_frac, t_ap_corr[TableColumn.EE_FRACTION],
                t_ap_corr[TableColumn.RADIUS]))

    def __init__(
            self, radius: float, sky_in: float, sky_out: float,
            verbose: int | bool = 0) -> None:
        """
        Aperture photometry called by starbug

        :param radius: Pixel radius of photometric aperture
        :type radius: float
        :param sky_in: Pixel radius of inner sky annuli
        :type sky_in: float
        :param sky_out: Pixel radius of outer sky annuli
        :type sky_out: float
        :param verbose: Set whether to print verbose output information
        :type verbose: bool
        """
        if sky_in < radius:
            warn("Sky annulus radii must be larger than aperture radius.\n")
            sky_in = radius + 1

        if sky_in >= sky_out:
            warn("Sky annulus outer radii must be larger than the inner.\n")
            sky_out = sky_in + 1

        self.radius: float = radius
        self.sky_in: float = sky_in
        self.sky_out: float = sky_out
        self.catalogue: Table = Table(None)
        self.verbose: int | bool = verbose

    def __call__(
            self, image: np.ndarray,
            detections: Table,
            error: np.ndarray | None = None,
            dq_flags: np.ndarray | None = None,
            ap_corr: float = 1.0,
            sig_sky: float = 3.0) -> Table:
        """
        Forced aperture photometry on a list of detections are an
         `astropy.table.Table` with columns x_centroid y_centroid or x_0 y_0
        This will add extra columns into this table ap_flux ap_sky

        :param detections: Table with source positions containing
                           x_centroid, y_centroid columns
        :type detections: astropy.table.Table
        :param image: 2D Image data to run photometry on
        :type image: numpy.ndarray
        :param error: 2D Image array containing photometric error per pixel
        :type error: numpy.ndarray
        :param dq_flags: 2D Image array containing JWST data quality flags per
                         pixel
        :type dq_flags: numpy.ndarray
        :param ap_corr: Aperture correction to be applied to the flux
        :type ap_corr: float
        :param sig_sky: Sigma threshold above the median to clip from sky
                        apertures
        :type sig_sky: float
        :return: Photometry catalogue
        :rtype: astropy.table.Table
        """
        return self._run(image, detections, error, dq_flags, ap_corr, sig_sky)

    def _run(self, image: np.ndarray,
             detections: Table,
             error: np.ndarray | None = None,
             dq_flags: np.ndarray | None = None,
             ap_corr: float = 1.0,
             sig_sky: float = 3.0) -> Table:
        """
        Forced aperture photometry on a list of detections are an
         `astropy.table.Table` with columns x_centroid y_centroid or x_0 y_0
        This will add extra columns into this table ap_flux ap_sky

        :param detections: Table with source positions containing
                           x_centroid, y_centroid columns
        :type detections: astropy.table.Table
        :param image: 2D Image data to run photometry on
        :type image: numpy.ndarray
        :param error: 2D Image array containing photometric error per pixel
        :type error: numpy.ndarray
        :param dq_flags: 2D Image array containing JWST data quality flags per
                         pixel
        :type dq_flags: numpy.ndarray
        :param ap_corr: Aperture correction to be applied to the flux
        :type ap_corr: float
        :param sig_sky: Sigma threshold above the median to clip from sky
                        apertures
        :type sig_sky: float
        :return: Photometry catalogue
        :rtype: astropy.table.Table
        """
        pos: list[tuple[Table | Row, Table | Row]]
        if (len({TableColumn.X_CENTROID, TableColumn.Y_CENTROID} &
                set(detections.colnames)) == 2):
            pos = [(line[TableColumn.X_CENTROID],
                    line[TableColumn.Y_CENTROID]) for line in detections]
        elif (len({TableColumn.X_0, TableColumn.Y_0} &
                  set(detections.colnames)) == 2):
            pos = [
                (line[TableColumn.X_0],
                 line[TableColumn.Y_0]) for line in detections]
        elif (len({TableColumn.X_INIT, TableColumn.Y_INIT} &
                  set(detections.colnames)) == 2):
            pos = [
                (line[TableColumn.X_INIT],
                 line[TableColumn.Y_INIT]) for line in detections]
        else:
            p_error(
                "Cannot identify position in detection catalogue ("
                "x_0/x_centroid)\n")
            raise Exception(
                "Cannot identify position in detection catalogue ("
                "x_0/x_centroid)\n")

        mask: np.ndarray = np.isnan(image)
        if error is None:
            error = np.sqrt(image)

        apertures: CircularAperture = CircularAperture(pos, self.radius)
        smooth_apertures: CircularAperture = CircularAperture(
            pos, min(1.5 * self.radius, self.sky_in))
        annulus_aperture: CircularAnnulus = CircularAnnulus(
            pos, r_in=self.sky_in, r_out=self.sky_out)

        self.log(
            "-> apertures: %.2g (%.2g - %.2g)\n" % (
                self.radius, self.sky_in, self.sky_out))
        phot: QTable = aperture_photometry(
            image, (apertures, smooth_apertures), error=error, mask=mask)
        self.catalogue = (
            Table(np.full((len(pos), 4), np.nan),
                  names=(TableColumn.SMOOTHNESS, TableColumn.FLUX,
                         TableColumn.E_FLUX, TableColumn.SKY)))

        self.log("-> calculating sky values\n")
        masks_raw = annulus_aperture.to_mask(method="center")

        # convert to list
        masks: list[ApertureMask] = (
            masks_raw if isinstance(masks_raw, list) else [masks_raw])

        # generate dat_list.
        dat_list: list[np.ndarray | None] = list(
            map(lambda a: a.multiply(image), masks))
        dat: np.ndarray

        try:
            dat = np.array(dat_list).astype(float)
        except (ValueError, TypeError) as e:
            # Cases where the array is inhomogeneous
            # If annulus reaches the edge of the image, it will create a
            # mask the wrong shape. If for whatever reason the point lies
            # outside the image, it will have None in the list, this needs
            # to be caught too
            warn(
                f"Ran into issues with the sky annuli. {e},"
                f" trying to fix them..\n")
            size: int = int(
                np.max([np.shape(d) for d in dat_list if d is not None]))
            fixed_dat: list[np.ndarray] = []
            for d in dat_list:
                if d is None:
                    fixed_dat.append(np.zeros((size, size)))
                elif (shape := np.shape(d)) != (size, size):
                    padded: np.ndarray = np.zeros((size, size))
                    padded[:shape[0], :shape[1]] += d
                    fixed_dat.append(padded)
                else:
                    fixed_dat.append(d)
            dat = np.array(fixed_dat)

        mask = (dat > 0 & np.isfinite(dat))
        dat[~mask] = np.nan
        clipped_dat: np.ma.MaskedArray = np.ma.MaskedArray(sigma_clip(
            dat.reshape(dat.shape[0], -1), sigma=sig_sky, axis=1))
        self.catalogue[TableColumn.SKY] = (
            np.ma.median(clipped_dat, axis=1).filled(fill_value=0))
        std: np.ndarray = np.ma.std(clipped_dat, axis=1)

        e_poisson: np.ndarray = phot[QTableColNames.SUM_ERR_0]
        esky_scatter: np.ndarray = apertures.area * std ** 2
        esky_mean: np.ndarray = (
            (std ** 2 * apertures.area ** 2) / annulus_aperture.area)

        self.catalogue[TableColumn.E_FLUX] = np.sqrt(
            e_poisson ** 2 + esky_scatter ** 2 + esky_mean ** 2)
        self.catalogue[TableColumn.FLUX] = (
            ap_corr * (phot[QTableColNames.SUM_0] - (
                self.catalogue[TableColumn.SKY] * apertures.area)))

        self.catalogue[TableColumn.FLUX][
            self.catalogue[TableColumn.FLUX] == 0] = np.nan

        ######################
        # Source "smoothness", the gradient of median pixel values within the
        # two test apertures
        ######################

        self.catalogue[TableColumn.SMOOTHNESS] = (
            (phot[QTableColNames.SUM_1] / smooth_apertures.area)
            / (phot[QTableColNames.SUM_0] / apertures.area))

        col: Column = Column(
            np.full(len(apertures), SourceFlags.SRC_GOOD),
            dtype=np.uint16, name=TableColumn.FLAG)
        if dq_flags is not None:
            self.log("-> flagging unlikely sources\n")
            for i, ap_mask in enumerate(apertures.to_mask(method="center")):
                i: int
                ap_mask: ApertureMask
                tmp: np.ndarray = ap_mask.multiply(dq_flags)
                if tmp is not None:
                    dq_dat = np.array(tmp, dtype=np.uint32)
                    if np.sum(dq_dat & (
                            DQFlags.DQ_DO_NOT_USE | DQFlags.DQ_SATURATED)):
                        col[i] |= SourceFlags.SRC_BAD
                    if np.sum(dq_dat & DQFlags.DQ_JUMP_DET):
                        col[i] |= SourceFlags.SRC_JMP
        self.catalogue.add_column(col)
        return self.catalogue

    def log(self, msg: str) -> None:
        """
        log message if in verbose mode

        :param msg: message to log
        :return: None
        """
        if self.verbose:
            printf(msg)
            sys.stdout.flush()
