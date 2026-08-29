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
from multiprocessing import shared_memory
from multiprocessing.shared_memory import SharedMemory

import numpy as np

from starbug2.command_line_interfaces.ast import ast_main
from starbug2.constants import ExitStates
from tests.generic import TEST_IMAGE_FITS, TEST_PSF_FITS


def run_ast_main(*args: str) -> int:
    """
    Helper function to wrap ast_main execution with predictable sequence
    arguments.

    :param args: Elements to append to the executable command array
    :return: The exit state integer from the binary
    """
    c: np.ndarray = np.array([0, 0, 0], dtype=np.int64)
    share_memory: SharedMemory = (
        shared_memory.SharedMemory(create=True, size=c.nbytes))
    loading_buffer: np.ndarray = np.ndarray(
        c.shape, dtype=c.dtype, buffer=share_memory.buf)
    return ast_main(
        ["starbug2-afs"] + list(args), share_memory, loading_buffer)


def test_ast_main_execution_states() -> None:
    """
    Verify that ast_main correctly handles valid image inputs and expected
    failures.
    """
    # Test behaviour with a valid, verified target image path
    assert (run_ast_main(
        TEST_IMAGE_FITS, "-s", "FILTER=F444W", "-s",
        f"PSF_FILE={TEST_PSF_FITS}") == ExitStates.EXIT_SUCCESS)

    # Test behaviour with an unrecognised, non-existent target path
    assert run_ast_main("nope") == ExitStates.EXIT_FAIL
