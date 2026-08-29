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

from typing import Final, Dict
from starbug2.constants import NIRCAM, DetectorLengths, STAR_BUG_MIRI


# (struct) containing JWST filter info
class FilterStruct:
    # noinspection SpellCheckingInspection, PyPep8Naming
    def __init__(
            self, full_width_half_max: float, instr: int, length: int,
            nlambda: int | None = None):
        """

        :param full_width_half_max:
        :param instr:
        :param length:
        :param nlambda: How many wavelengths to model for broadband?
            The default depends on how wide the filter is: (5,3,1) for
            types (W,M,N) respectively in webb.stpsf, and its default is 40
        """
        self._full_width_half_max: float = full_width_half_max
        self._instr: int = instr
        self._length: int = length
        self._nlambda: int | None = nlambda

    @property
    def full_width_half_max(self) -> float:
        return self._full_width_half_max

    @property
    def instr(self) -> int:
        return self._instr

    @property
    def length(self) -> int:
        return self._length

    # noinspection SpellCheckingInspection, PyPep8Naming
    @property
    def nlambda(self) -> int | None:
        return self._nlambda


# as of 08/06/2023
STAR_BUG_FILTERS: Final[Dict[str, FilterStruct]] = {
    "F070W": FilterStruct(0.742, NIRCAM, DetectorLengths.SHORT),
    "F090W": FilterStruct(0.968, NIRCAM, DetectorLengths.SHORT),
    "F115W": FilterStruct(1.194, NIRCAM, DetectorLengths.SHORT),
    "F140M": FilterStruct(1.484, NIRCAM, DetectorLengths.SHORT),
    "F150W": FilterStruct(1.581, NIRCAM, DetectorLengths.SHORT),
    "F162M": FilterStruct(1.710, NIRCAM, DetectorLengths.SHORT),
    "F164N": FilterStruct(1.742, NIRCAM, DetectorLengths.SHORT),
    # NOTE: MAJOR CHANGE
    # the 20 here is to resolve an issue inside stpsf calc_psf for F150W2 as
    # the wave max value being  2.35 * 1e-6 and the wave lengths when
    # partitioned at 40 results in a final wave length of
    # 2.3626144323647297e-06 and fails validation. By making nlambda 20, ive
    # added more entries into each bucket, which I think just reduces the
    # samples, but keeps the structure of the PSF the same
    "F150W2": FilterStruct(1.452, NIRCAM, DetectorLengths.SHORT, 20),
    "F182M": FilterStruct(1.935, NIRCAM, DetectorLengths.SHORT),
    "F187N": FilterStruct(1.968, NIRCAM, DetectorLengths.SHORT),
    "F200W": FilterStruct(2.065, NIRCAM, DetectorLengths.SHORT),
    "F210M": FilterStruct(2.194, NIRCAM, DetectorLengths.SHORT),
    "F212N": FilterStruct(2.226, NIRCAM, DetectorLengths.SHORT),
    "F250M": FilterStruct(1.302, NIRCAM, DetectorLengths.LONG),
    "F277W": FilterStruct(1.397, NIRCAM, DetectorLengths.LONG),
    "F300M": FilterStruct(1.540, NIRCAM, DetectorLengths.LONG),
    "F322W2": FilterStruct(1.524, NIRCAM, DetectorLengths.LONG),
    "F323N": FilterStruct(1.683, NIRCAM, DetectorLengths.LONG),
    "F335M": FilterStruct(1.730, NIRCAM, DetectorLengths.LONG),
    "F356W": FilterStruct(1.810, NIRCAM, DetectorLengths.LONG),
    "F360M": FilterStruct(1.873, NIRCAM, DetectorLengths.LONG),
    "F405N": FilterStruct(2.095, NIRCAM, DetectorLengths.LONG),
    "F410M": FilterStruct(2.111, NIRCAM, DetectorLengths.LONG),
    "F430M": FilterStruct(2.206, NIRCAM, DetectorLengths.LONG),
    "F444W": FilterStruct(2.222, NIRCAM, DetectorLengths.LONG),
    "F460M": FilterStruct(2.397, NIRCAM, DetectorLengths.LONG),
    "F466N": FilterStruct(2.413, NIRCAM, DetectorLengths.LONG),
    "F470N": FilterStruct(2.444, NIRCAM, DetectorLengths.LONG),
    "F480M": FilterStruct(2.492, NIRCAM, DetectorLengths.LONG),
    "F560W": FilterStruct(1.882, STAR_BUG_MIRI, DetectorLengths.NULL),
    "F770W": FilterStruct(2.445, STAR_BUG_MIRI, DetectorLengths.NULL),
    "F1000W": FilterStruct(2.982, STAR_BUG_MIRI, DetectorLengths.NULL),
    "F1130W": FilterStruct(3.409, STAR_BUG_MIRI, DetectorLengths.NULL),
    "F1280W": FilterStruct(3.818, STAR_BUG_MIRI, DetectorLengths.NULL),
    "F1500W": FilterStruct(4.436, STAR_BUG_MIRI, DetectorLengths.NULL),
    "F1800W": FilterStruct(5.373, STAR_BUG_MIRI, DetectorLengths.NULL),
    "F2100W": FilterStruct(6.127, STAR_BUG_MIRI, DetectorLengths.NULL),
    "F2550W": FilterStruct(7.300, STAR_BUG_MIRI, DetectorLengths.NULL),
}
