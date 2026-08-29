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
from typing import List, Any

import numpy as np
from astropy.io.fits import PrimaryHDU, ImageHDU, BinTableHDU
from astropy.visualization import ZScaleInterval
from astropy.table import Row, Table
from scipy.interpolate import RegularGridInterpolator
import matplotlib.image as mpimg
from matplotlib.colors import LinearSegmentedColormap
from astropy.wcs import WCS
from starbug2.constants import URL_DOCS, HeaderTags, TableColumn
from starbug2.utilities import utils
from starbug2.utilities.filters import STAR_BUG_FILTERS

# try to import pyplot as plt.
try:
    import matplotlib.pyplot as plt
except ImportError:
    from matplotlib import use
    use("TkAgg")
    import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure


def load_style(f_name: str) -> None:
    """
    Load a pyplot style sheet

    :param f_name: Filename of the style sheet
    :type f_name: str
    :return: None
    """
    if os.path.exists(f_name):
        plt.style.use(f_name)
    else:
        utils.p_error("Unable to load style sheet \"%s\"\n" % f_name)


def _generate_regular_grid_interpolator(
        x, y, bins) -> RegularGridInterpolator:
    """
    generates a regular grid interpolator

    :param x: x coord
    :param y: y coord
    :param bins: the bins.
    :return: the configured RegularGridInterpolator
    :rtype: a RegularGridInterpolator
    """
    hist: np.ndarray
    xx: np.ndarray
    yy: np.ndarray
    hist, _, _ = np.histogram2d(x, y, bins=bins)
    xx = np.linspace(min(x), max(x), bins)
    yy = np.linspace(min(y), max(y), bins)
    return RegularGridInterpolator((xx, yy), hist)


def plot_test(axes: Axes) -> Axes:
    """
    Just plot the starbug image

    :param axes: Axis to plot into
    :type axes: plt.axes
    :return: The working axes
    :rtype: plt.axes
    """
    if not utils.wget(URL_DOCS, f_name="/tmp/sb.png"):
        axes.imshow(mpimg.imread("/tmp/sb.png"))
        axes.set_title("starbug2 v%s" % utils.get_version())
        axes.set_axis_off()
    return axes


def plot_cmd(
        tab: Table,
        colour: str,
        mag: str,
        axis: Axes | None = None,
        col: str | tuple[float, ...] | None = None,
        hess: bool = True,
        x_lim: tuple[float, float] | None = None,
        y_lim: tuple[float, float] | None = None,
        **kwargs: Any) -> Axes:
    """
    Plot a Colour-Magnitude Diagram (CMD) with an optional Hess-based density
    colouring.

    :param tab: Astropy Table containing the stellar catalogue data
    :type tab: astropy.table.Table
    :param colour: The column name representing the colour index
                   (e.g., 'B-V' or 'g-r')
    :type colour: str
    :param mag: The column name representing the magnitude (e.g., 'V' or 'g')
    :type mag: str
    :param axis: The matplotlib Axes object to plot on, defaults to None
    :type axis: matplotlib.axes.Axes or None
    :param col: Custom colour string or sequence for the scatter
                points/colormap gradient
    :type col: str or tuple or None
    :param hess: Whether to calculate and apply a density-based Hess plot
                 colour gradient
    :type hess: bool
    :param x_lim: X-axis bounds for the colour index (min, max)
    :type x_lim: tuple of (float, float) or None
    :param y_lim: Y-axis bounds for the magnitude (min, max)
    :type y_lim: tuple of (float, float) or None
    :param kwargs: Additional keyword arguments passed directly to ax.scatter
    :type kwargs: dict
    :return: The matplotlib Axes object with the rendered plot
    :rtype: matplotlib.axes.Axes
    """
    tt: Table = utils.colour_index(tab, [colour, mag])
    mask: np.ndarray = ~ (tt[colour].mask | tt[mag].mask)
    cc: np.ndarray = tt[colour][mask]
    mm: np.ndarray = tt[mag][mask]

    if not axis:
        _, axis = plt.subplots(1)

    if x_lim is None:
        x_lim = (float(np.nanmin(cc)), float(np.nanmax(cc)))
    if y_lim is None:
        y_lim = (np.nanmin(mm), np.nanmax(mm))

    spatial_mask: np.ndarray = (
        (cc >= x_lim[0]) & (cc <= x_lim[1]) &
        (mm >= y_lim[0]) & (mm <= y_lim[1]))
    cc = cc[spatial_mask]
    mm = mm[spatial_mask]

    # apply default colour
    if col is None:
        col = plt.rcParams["axes.prop_cycle"].by_key()["color"][0]

    # make segmented colour map
    cmap: LinearSegmentedColormap = LinearSegmentedColormap.from_list(
        "", [plt.rcParams["axes.prop_cycle"].by_key()["color"][0], col])

    if hess:
        bins: int = 100
        f: RegularGridInterpolator = (
            _generate_regular_grid_interpolator(cc, mm, bins))
        col = [f([X, Y]) for X, Y in zip(cc, mm)]
    pyplot_kw: dict[str, int] = {"lw": 0, "s": 3}
    pyplot_kw.update(kwargs)
    axis.scatter(cc, mm, c=col, cmap=cmap, **pyplot_kw)

    axis.set_xlabel(colour)
    axis.set_ylabel(mag)
    axis.set_xlim(x_lim)

    # Invert the Y-axis because brighter astronomical magnitudes have
    # lower values
    axis.set_ylim(*y_lim[::-1])
    return axis


def plot_inspect_source(
        src: Row, images: List[PrimaryHDU | ImageHDU | BinTableHDU | None]):
    """
    Show a source in an array of images

    :param src: Input source to look at
    :type src: astropy.Table Row
    :param images: List of fits images to inspect
    :type images: HDUList
    :return: The figure
    :rtype: plt.figure
    """

    n: int = len(images)
    figure: Figure
    axs: Axes | List[Axes]
    figure, axs = plt.subplots(1, n, figsize=(1.7 * n, 2))
    if n == 1:
        axs = [axs]  # noqa
    images: List[ImageHDU | PrimaryHDU | BinTableHDU | None] = sorted(
        images, key=lambda a:
            list(STAR_BUG_FILTERS.keys()).index(a.header[HeaderTags.FILTER]))

    # arcsec?
    size: float = 0.1
    n: int
    im: ImageHDU | PrimaryHDU
    axis: Axes
    assert isinstance(axs, list)
    for n, (im, axis) in enumerate(zip(images, axs)):
        wcs: WCS = WCS(im)
        x: np.ndarray
        y: np.ndarray
        x, y = wcs.all_world2pix(
            np.ones(2) * src["RA"],
            np.array([1, 1 + (size / 3600)]) * src["DEC"], 0)
        dp: np.ndarray = np.sqrt((x[1] - x[0]) ** 2 + (y[1] - y[0]) ** 2)

        x_min: int = max(0, int(np.floor(x[0] - dp)))
        x_max: int = min(im.data.shape[1] - 1, int(np.ceil(x[0] + dp)))
        y_min: int = max(0, int(np.floor(y[0] - dp)))
        y_max: int = min(im.data.shape[0] - 1, int(np.ceil(y[0] + dp)))

        dat: np.ndarray = im.data[
              min(y_min, y_max): max(y_min, y_max),
              min(x_min, x_max): max(x_min, x_max)]
        if all(dat.shape):
            axis.imshow(ZScaleInterval()(dat), cmap="Greys_r", origin="lower")
            axis.text(0, 0, im.header.get(HeaderTags.FILTER), c="white")

        axis.set_axis_off()
        figure.suptitle(src[TableColumn.CAT_NUM][0])
    figure.tight_layout()
    return figure
