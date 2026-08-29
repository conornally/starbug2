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
import numpy as np
from typing import cast, Any, Tuple, Callable, Dict

from astropy.io.fits import ImageHDU, PrimaryHDU
from photutils.datasets import make_model_image, make_random_models_table
from astropy.table import Table, QTable, Column
from astropy.io import fits
from photutils.psf import ImagePSF
from scipy.optimize import curve_fit
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from starbug2.core.star_bug_config import StarBugMainConfig
from starbug2.constants import ExitStates, TableColumn, ERR

try:
    import matplotlib.pyplot as plt
except ImportError:
    import matplotlib
    matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt

from starbug2.utilities.utils import (
    printf, p_error, get_mj_ysr2jy_scale_factor, warn, flux_to_pogson_mag,
    ext_names)


class ArtificialStars:

    @staticmethod
    def verify(
            config: StarBugMainConfig,
            main_image_shape: Tuple[int, int]) -> ExitStates:
        if (config.test_magnitude_bright_limit -
                config.test_magnitude_faint_limit >= 0):
            warn("Detected magnitude range in wrong order,"
                 " put bright limit first\n")
            return ExitStates.EXIT_FAIL

        base_shape: np.ndarray = np.copy(main_image_shape)
        if any(base_shape < config.sub_image_crop_size):
            p_error("sub image size greater than image size, setting to "
                    "'safe' value %d.\n" % config.sub_image_crop_size)
        return ExitStates.EXIT_SUCCESS

    @staticmethod
    def add_stars(
            base_image: fits.HDUList | None, config: StarBugMainConfig,
            buffer: int, psf: np.ndarray | None, n_hdu: int) -> Tuple[
                ExitStates, QTable, fits.HDUList]:
        """
        adds new stars to the image.
        :param base_image: copy of the current image
        :type base_image: fits.HDUList or None
        :param config: the main config
        :type config: StarBugMainConfig
        :param buffer: the buffer
        :type buffer: int
        :param psf: the point source function
        :type psf: np.ndarray | None
        :param n_hdu: the n_hdu.
        :type n_hdu: int
        :return: exit state success the location of the new stars and the
                 new image
        :rtype: Tuple[ExitStates, atrophy.QTable, fits.HDUList]
        """

        # this line is due to the fits file being a lazy reader. so this is not
        # in memory, it is still accessing the file directly. So a copy avoids
        # corrupting the original file.
        assert base_image is not None
        image: fits.HDUList = base_image.__deepcopy__()

        shape: tuple[int, int] = image[n_hdu].shape # noqa
        main_image: ImageHDU | PrimaryHDU = image[n_hdu]

        # create list of fake star locations and add new magnitudes
        source_list: QTable = make_random_models_table(
            config.stars_per_artificial_test, {
                TableColumn.X_0: [buffer, shape[0] - buffer],
                TableColumn.Y_0: [buffer, shape[1] - buffer],
                TableColumn.MAG:
                    [config.test_magnitude_bright_limit,
                     config.test_magnitude_faint_limit]
            }, config.ast_seed
        )

        # utilising the Pogsons brightness equation. it ensures an astronomical
        # magnitude is converted to a linear flux density if the zero point
        # magnitude is correct in that it "defined as the magnitude of a
        # source that produces exactly $1\text{ Jy}$ of flux " see
        # ref https://adsabs.harvard.edu/pdf/1983ApJ...266..713O#
        source_list.add_column(
            10.0 ** (
                (source_list[TableColumn.MAG] - config.zero_point_magnitude)
                / -2.5),
            name=TableColumn.FLUX)
        source_list.remove_column(TableColumn.ID)

        # create new image bits for each fake star.
        scale_factor: float | int = (
            get_mj_ysr2jy_scale_factor(main_image))
        image_psf = ImagePSF(psf)

        # for sanity sakes. let's run them 1 at a time and get the difference
        # between what flux is asked for and what's added.
        retrieved_flux: list = []
        retrieved_mag: list = []
        retrieved_mag_diff: list = []
        for i in range(len(source_list)):
            single_source_table = source_list[i: i + 1]

            star_overlay: np.ndarray = (
                make_model_image(
                    shape, image_psf, single_source_table,
                    model_shape=image_psf.data.shape)
                / scale_factor)

            # apply the new data to the original image.
            main_image_new_data = main_image.data + star_overlay
            image[n_hdu].data = main_image_new_data

            # add extra columns
            total_star_flux = float(np.sum(star_overlay * scale_factor))
            retrieved_flux.append(total_star_flux)
            added_mag, _ = flux_to_pogson_mag(total_star_flux)
            added_mag = added_mag + config.zero_point_magnitude
            retrieved_mag.append(added_mag)
            retrieved_mag_diff.append(
                float(single_source_table[TableColumn.MAG][0]) - added_mag)

        if ERR in ext_names(image):
            image[ERR].data += np.sqrt(np.abs(main_image.data))

        source_list[TableColumn.TOTAL_FLUX_ADDED] = retrieved_flux
        source_list[TableColumn.TOTAL_MAG] = retrieved_mag
        source_list[TableColumn.TOTAL_DIFF_MAG] = retrieved_mag_diff

        # if told to generate the added list. generate.
        if config.ast_save_added_image:
            os.makedirs(config.ast_save_added_image_path, exist_ok=True)
            image.verify('fix')
            image.writeto(os.path.join(
                config.ast_save_added_image_path,
                f"inserted_image_for_test_{config.ast_test_index}.fits"))
            source_list.write(os.path.join(
                config.ast_save_added_image_path,
                f"stars_data_for_test_{config.ast_test_index}.fits"))

        # the rest of the code don't want these values, so delete them from the
        # source list.
        source_list.remove_column(TableColumn.TOTAL_FLUX_ADDED)
        source_list.remove_column(TableColumn.TOTAL_MAG)
        source_list.remove_column(TableColumn.TOTAL_DIFF_MAG)

        # hand back the source list, the new image, and the exit state.
        return ExitStates.EXIT_SUCCESS, source_list, image


def get_completeness(test_result: Table) -> Table:
    """
    Compile the results into magnitude binned values of recovery fraction
    and flux error.

    :param test_result: The output from auto_run.
    :type test_result: astropy.table.Table
    :return: A table containing percentage completeness as a function of
             magnitude.
    :rtype: astropy.table.Table
    """

    bins: np.ndarray = np.arange(
        np.floor(np.nanmin(test_result[TableColumn.MAG])),
        np.ceil(np.nanmax(test_result[TableColumn.MAG])),
        0.1)
    percents: np.ndarray = np.zeros(len(bins))
    errors: np.ndarray = np.zeros(len(bins))
    offsets: np.ndarray = np.zeros(len(bins))
    means: np.ndarray = np.zeros(len(bins))

    i_bins: np.ndarray = np.atleast_1d(np.digitize(
        test_result[TableColumn.MAG], bins=bins))
    for i in range(1, int(max(bins)) + 1):
        indices = np.where(i_bins == i)[0]
        # skip if no data.
        if len(indices) == 0:
            continue

        # process data
        binned: Table = test_result[indices]
        if len(binned) > 0:
            percents[i] = float(
                sum(binned[TableColumn.STATUS])) / len(binned)

        mag_inj: np.ndarray
        mag_det: np.ndarray

        # as were comparing differences. we do not need to worry about
        # Zero Point.
        mag_inj, _ = flux_to_pogson_mag(binned[TableColumn.FLUX])
        mag_det, _ = flux_to_pogson_mag(binned[TableColumn.FLUX_DET])

        mag_difference: np.ndarray = mag_inj - mag_det
        errors[i] = np.nanstd(mag_difference)
        means[i] = np.nanmean(mag_difference)
        offsets[i] = np.nanmedian(
            binned[TableColumn.FLUX] / binned[TableColumn.FLUX_DET])

    out: Table = Table(
        [bins, percents, errors, offsets],
        names=(TableColumn.MAG, TableColumn.REC, TableColumn.ERR_LOWER,
               TableColumn.OFF),
        dtype=(float, float, float, float))
    return out


def get_spatial_completeness(
        test_result: Table, image: np.ndarray | None,
        res: int = 10) -> np.ndarray | None:
    """
    Produce an image array showing the spatially dependent recovery fraction.

    :param test_result: The output from auto_run.
    :type test_result: astropy.table.Table
    :param image: 2D image array to take the shape from.
    :type image: numpy.ndarray
    :param res: The resolution of the spatial bins.
    :type res: int
    :return: A 2D array the same shape as the image input, where pixel values
        show the fraction of injected sources recovered in this bin.
    :rtype: numpy.ndarray
    """
    if image is None:
        return None

    x_bins: np.ndarray = np.arange(
        min(np.atleast_1d(test_result[TableColumn.X_0])),
        max(np.atleast_1d(test_result[TableColumn.X_0])), int(res))
    y_bins: np.ndarray = np.arange(
        min(np.atleast_1d(test_result[TableColumn.Y_0])),
        max(np.atleast_1d(test_result[TableColumn.Y_0])), int(res))
    percents: np.ndarray = np.zeros(image.shape)

    xi: int
    for xi in x_bins[:-1]:
        xo: int = xi + res
        yi: int
        for yi in y_bins[:-1]:
            yo: int = yi + res
            mask: np.ndarray = (
                (test_result[TableColumn.X_0] >= xi) &
                (test_result[TableColumn.X_0] < xo) &
                (test_result[TableColumn.Y_0] >= yi) &
                (test_result[TableColumn.Y_0] < yo))
            binned: Table = test_result[mask]
            if len(binned):
                percents[int(yi): int(yo), int(xi): int(xo)] = (
                    float(np.sum(binned[TableColumn.STATUS])) / len(binned))
    return percents


def estimate_completeness_mag(ast: Table) -> (
        Tuple[Tuple[float, float, float] | None,
              Tuple[float, float, float] | None]):
    """
    Estimate the completeness level of the artificial star test.

    :param ast: Output of Artificial_Stars.get_completeness, table must
               contain columns (mag, rec).
    :type ast: astropy.table.Table
    :return: A tuple containing:
        - **fit** (*list*): The fitting parameters to the logistic curve
                            $f(x) = \frac{l}{1 + exp(-k(x - x_0))}$ formatted
                            as ``[l, x_0, k]``.
        - **complete** (*list*): Magnitude of 90%, 70% and 50% completeness.
    :rtype: tuple[list, list]
    """
    fit: Tuple[float, float, float] | None = None
    completeness: Tuple[float, float, float] | None = None

    # Syntax: Callable[[Param1Type, Param2Type, ...], ReturnType]
    fn_i: Callable[[float, float, float, float], float] = (
        lambda y, limit, k, xo: xo - (np.log((limit / y) - 1) / k)
    )

    if len(set(ast.colnames) & {TableColumn.MAG, TableColumn.REC}) == 2:
        try:
            # need the *_ as the return tuple can be multiple sizes. The *_
            # allows the IDE to not freak out, especially as we don't care
            # about the rest of the return values.
            bounds = ([0.8, -np.inf, 0], [1.0, np.inf, np.inf])
            fit, *_ = curve_fit(
                scurve, ast[TableColumn.MAG], ast[TableColumn.REC],
                [1, -1, np.median(ast[TableColumn.MAG])],
                bounds=bounds)
            assert fit is not None
            completeness = (fn_i(0.9, *fit), fn_i(0.7, *fit), fn_i(0.5, *fit))
        except (RuntimeError, ValueError) as e:
            warn(f"Unable to fit completeness fractions: {e}\n")
    else:
        p_error("Input table must have columns 'mag' and 'rec'\n")
    return fit, completeness


def scurve(
        x: np.ndarray, limit: float, k: float,
        xo: float) -> float | np.ndarray:
    """
    S-curve function to fit completeness results to.

    math::  f(x) = \\frac{l}(1 + exp(-k(x - x_0)))

    :param x: Magnitude range or array to input into the function.
    :type x: list or numpy.ndarray
    :param limit: Maximum value asymptote (typically representing maximum
              completeness, near 1.0).
    :type limit: float
    :param xo: The inflection point of the curve (the magnitude where
               completeness is 50%).
    :type xo: float
    :param k: The logistic growth rate or steepness of the curve.
    :type k: float
    :return: Calculated function value(s) matching the shape of the
             input ``x``.
    :rtype: float or numpy.ndarray
    """
    return limit / (1 + np.exp(-k * (x - xo)))


def plot_mid_plot(
        ax: Axes, mag_true: Column, delta_mag, bin_centers, med_offsets,
        std_offsets) -> None:
    ax.scatter(
        mag_true,
        delta_mag,
        alpha=0.05,
        c='gray',
        s=2,
        label='Individual ASTs'
    )

    # Reference line at 0 offset (no bias)
    ax.axhline(0, color='black', linestyle='--', linewidth=1, alpha=0.7)

    # Binned median trend line
    ax.plot(
        bin_centers,
        med_offsets,
        color='red',
        linewidth=2,
        label='Median Offset'
    )

    # Error band representing spread
    ax.fill_between(
        bin_centers,
        med_offsets - std_offsets,
        med_offsets + std_offsets,
        color='red',
        alpha=0.2
    )

    # noinspection SpellCheckingInspection
    ax.set_ylabel(r'$\Delta m$ ($m_{\mathrm{det}} - m_{\mathrm{inj}}$)')
    # Zoom in on the bias zone (-0.5 to +0.5 mag is standard)
    ax.set_ylim(-0.5, 0.5)
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(loc='upper left')


def plot_top_plot(
        ax: Axes, completeness_raw: Table,
        completeness: Tuple[float, float, float],
        cfit: Tuple[float, float, float], filter_string: str | None,
        plot_ast: str) -> None:
    ax.scatter(
        completeness_raw[TableColumn.MAG],
        completeness_raw[TableColumn.REC], c='k', lw=0, s=8)
    ax.plot(completeness_raw[TableColumn.MAG],
            scurve(completeness_raw[TableColumn.MAG], *cfit),
            c='g',
            label=r"$f(x)=\frac{%.2f}{1+e^{%.2f("r"x-%.2f)}}$" % (
                cfit[0], cfit[1], cfit[2]))
    ax.axvline(
        completeness[0], c="seagreen", ls='--',
        label=("90%%:%.2f" % completeness[0]), lw=0.75)
    ax.axvline(
        completeness[1], c="seagreen", ls='-.',
        label=("70%%:%.2f" % completeness[1]), lw=0.75)
    ax.axvline(
        completeness[2], c="seagreen", ls=':',
        label=("50%%:%.2f" % completeness[2]), lw=0.75)
    ax.scatter(completeness, (0.9, 0.7, 0.5), marker='*', c='teal', s=10)
    ax.tick_params(direction="in", top=True, right=True)
    ax.set_title("Artificial Star Test")

    assert filter_string is not None
    ax.set_xlabel(filter_string)
    ax.set_ylabel("Fraction Recovered")
    ax.set_yticks([0, .25, .5, .75, 1])
    ax.legend(loc="lower left", frameon=False, fontsize=8)
    plt.tight_layout()
    printf("--> %s\n" % plot_ast)


def photometric_bias(
        raw: Table) -> Tuple[Table, np.ndarray, np.ndarray, np.ndarray]:
    """
    generate the photometric bias data
    :param raw: the result table.
    :type raw: table.
    :return: the differences in magnitude, the centres of each bin,
    the medium offsets, the standard deviation offsets.
    :rtype: Tuple[Table, np.ndarray, np.ndarray, np.ndarray]
    """
    bins: np.ndarray = np.arange(
        np.floor(np.nanmin(raw[TableColumn.MAG])),
        np.ceil(np.nanmax(raw[TableColumn.MAG])),
        0.1)
    mag_det, _ = flux_to_pogson_mag(raw[TableColumn.FLUX_DET])

    delta_mag: Table = raw[TableColumn.MAG] - mag_det
    med_offsets = np.full(len(bins) - 1, np.nan)
    std_offsets = np.full(len(bins) - 1, np.nan)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    for i in range(1, len(bins)):
        # Select stars in the current magnitude bin
        mask = ((raw[TableColumn.MAG] >= bins[i - 1]) &
                (raw[TableColumn.MAG] < bins[i]))
        bin_deltas = delta_mag[mask]
        valid_deltas = bin_deltas[np.isfinite(bin_deltas)]

        if len(valid_deltas) > 0:
            med_offsets[i - 1] = np.nanmedian(valid_deltas)
        # Use 16th to 84th percentile / 2 or std for error bounds
        std_offsets[i - 1] = np.nanstd(valid_deltas)
    return delta_mag, bin_centers, med_offsets, std_offsets


def _generate_head(
        completeness: Tuple[float, float, float],
        cfit: Tuple[float, float, float],
        filter_string: str | None) -> Dict[str, str | float]:
    """
    generates the head of the results table.
    :param completeness: the completeness
    :type completeness:  Tuple[float, float, float]
    :param cfit: the coefficients.
    :type cfit:  Tuple[float, float, float]
    :param filter_string:  the filter string.
    :type filter_string: str | None
    :return:
    """
    head: Dict[str, str | float] = {
        "COMPLETE_FN": "F(x)=l/(1+exp(-k(x-xo)))", "l": cfit[0],
        "k": cfit[1], "xo": cfit[2]}
    for i, frac in enumerate((90, 70, 50)):
        if completeness[i] and not np.isnan(completeness[i]):
            printf(
                "-> complete to %d%%: %s=%.2f\n" % (
                    frac, filter_string, completeness[i]))
            head["COMPLETE %d%%" % frac] = str(completeness[i])
    return head


def add_mag_columns(raw: Table, config: StarBugMainConfig) -> Table:
    """
    adds extra columns for astrophysics's to use so they could generate
     their own plots.
    :param raw: the raw table.
    :param config: the config object.
    :type config: StarBugMainConfig
    :return: a new table with MAG detection and MAG difference columns
             inserted.
    :rtype: Table
    """
    mag_det: np.ndarray
    mag_det, _ = flux_to_pogson_mag(raw[TableColumn.FLUX_DET])
    mag_det += config.zero_point_magnitude
    mag_difference: np.ndarray
    mag_difference: np.ndarray = raw[TableColumn.MAG] - mag_det

    raw_data: Table = raw.copy()
    raw_data.add_column(mag_det, name=TableColumn.MAG_DET)
    raw_data.add_column(mag_difference, name=TableColumn.MAG_DIFF)
    return raw_data


def compile_results(
        raw: Table, config: StarBugMainConfig,
        image: np.ndarray | None = None,
        plot_ast: str | None = None,
        filter_string: str | None = "m") -> dict[str, fits.HDUList]:
    """
    Compile all the raw data into usable results

    :param raw:raw data
    :type raw: astro.table.table
    :param image: the image data
    :type image: np.ndarray
    :param plot_ast: the save plot file name
    :type plot_ast: str or None
    :param filter_string: the filter string
    :type filter_string: str | None
    :param config: the config object
    :type config: StarBugMainConfig
    :return: the results
    :rtype: fits.HDUList
    """

    completeness_raw: Table = get_completeness(raw)
    cfit: Tuple[float, float, float]
    completeness: Tuple[float, float, float]
    cfit, completeness = estimate_completeness_mag(completeness_raw)
    spatial_completeness: np.ndarray | None = (
        get_spatial_completeness(raw, image, res=10))
    delta_mag, bin_centers, med_offsets, std_offsets = photometric_bias(raw)

    head: Dict[str, str | float]
    head = _generate_head(completeness, cfit, filter_string)

    # add mag out and mag difference columns to the raw table.
    mag_raw = add_mag_columns(raw, config)

    results: dict[str, fits.HDUList] = {
        "-ast.fits": fits.HDUList(
            [fits.PrimaryHDU(header=fits.Header(head)),
             fits.BinTableHDU(data=completeness_raw, name="AST"),
             fits.BinTableHDU(data=mag_raw, name="RAW")]),
        "-ast-spatial.fits": fits.HDUList(
            [fits.PrimaryHDU(header=fits.Header(head)),
             fits.ImageHDU(data=cast(Any, spatial_completeness), name="CMP")])}

    if plot_ast:
        fig: Figure
        ax1: Axes
        ax2: Axes
        ax3: Axes
        fig, (ax1, ax2, ax3) = plt.subplots(
            3, figsize=(8, 10), dpi=300, sharex=True)

        plot_top_plot(
            ax1, completeness_raw, completeness, cfit, filter_string,
            plot_ast)
        plot_mid_plot(
            ax2, raw[TableColumn.MAG], delta_mag, bin_centers, med_offsets,
            std_offsets)
        fig.savefig(plot_ast, dpi=300)
    return results
