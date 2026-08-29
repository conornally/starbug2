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
from multiprocessing.shared_memory import SharedMemory
from multiprocessing import Process, shared_memory
from typing import Tuple

import numpy as np
from time import sleep
from astropy.table import Table

from starbug2.core.main_components.multi_treading_execution import (
    execute_one_core_run_ast, execute_multicore_ast)
from starbug2.core.main_components.one_time_runs import ast_one_time_runs
from starbug2.constants import ExitStates
from starbug2.core.star_bug_config import StarBugMainConfig
from starbug2.core.starbug_main import StarbugBase
from starbug2.utilities.utils import printf, p_error, parse_cmd

import photutils

# Force photutils to strictly return standard QTables globally
photutils.future_column_names = True


def load(loading_buffer: np.ndarray) -> None:
    """
    A loading bar that should be run in a subprocess
    It sits and watches the shared memory buffer and periodically
    prints out a progress bar
    """
    while loading_buffer[0] < loading_buffer[1]:
        sleep(1)
        p: np.ndarray = loading_buffer[0] / loading_buffer[1]
        msg: str = f"recovering:{loading_buffer[2]}%"
        s: str = "\x1b[2K%s|%-40s|%d/%d\r" % (
            msg, int(p*40)*'=', int(loading_buffer[0]), int(loading_buffer[1]))
        printf(s)
        sys.stdout.flush()
    printf("\n")


def ast_parse_argv(argv: list[str]) -> StarBugMainConfig:
    """
    Organise the sys argv line into options, values and arguments

    :param argv: the arguments
    :return: the config class
    :rtype: StarBugMainConfig
    """
    config: StarBugMainConfig = StarBugMainConfig()
    short_definition: str
    long_definition: list[str]

    short_definition, long_definition = (
        config.generate_ast_get_opt_definitions())
    _, argv = parse_cmd(argv)
    config.populate_params(
        argv, short_definition, long_definition, config.AST_FLAG_MAP)
    return config


def report_before_logs(config: StarBugMainConfig, f_name: str) -> None:
    """
    does some verbose logging about what's going to happen.
    :param config: the starbug config
    :type config: StarBugMainConfig
    :param f_name: the first fits image file name
    :type f_name: str
    :return: None
    :rtype: None
    """
    printf("Artificial Stars\n----------------\n")
    printf("-> loading %s\n" % f_name)
    if config.param_file:
        printf("-> parameters: %s\n" % config.param_file)
    printf("-> running %d tests with %d injections per test\n" % (
        config.artificial_star_tests_count, config.stars_per_artificial_test))
    printf("-> magnitude range: %.1f - %.1f\n" % (
        config.test_magnitude_bright_limit,
        config.test_magnitude_faint_limit))
    if config.ast_no_psf_phot:
        printf("-> skipping PSF photometry step\n")
    if config.ast_no_background:
        printf("-> skipping background estimation step\n")


def start_buffer(loading_buffer: np.ndarray, n_tests: int) -> Process:
    """
    initializes the loading buffer.
    :param loading_buffer: the loading buffer.
    :type loading_buffer: np.ndarray.
    :param n_tests: the number of tests.
    :type n_tests: int
    :return: the loading process.
    :rtype: Process
    """
    loading_buffer[0] = 0
    loading_buffer[1] = n_tests
    loading: Process = Process(target=load, args=[loading_buffer])
    loading.start()
    return loading


def execute_ast(
        config: StarBugMainConfig, f_name: str,
        loading_buffer: np.ndarray) -> Tuple[list[Table | None], ExitStates]:
    """
    determine execution state in terms of processors available and execute.
    :param config: the main config
    :type config: StarBugMainConfig
    :param f_name: the first fits name.
    :type f_name: str
    :param loading_buffer: the loading buffer.
    :type loading_buffer: np.ndarray
    :return: the results
    :rtype: Tuple[list[Table | None], ExitStates]
    """
    # ensure we read in the psf and image objects and fill the config
    # accordingly
    psf_image_config: StarBugMainConfig = StarBugMainConfig()
    psf_image_config.ast_load_psf = True
    psf_image_config.custom_filter = config.custom_filter
    star_bug = StarbugBase(f_name, psf_image_config, None, None)
    star_bug.run_starbug()
    config.ast_psf = star_bug.psf

    # management
    if config.verbose_logs:
        report_before_logs(config, f_name)
    loading: Process = start_buffer(
        loading_buffer, config.artificial_star_tests_count)

    # Initialise output container tracking tables
    outs: list[Table | None]
    exit_state: ExitStates
    results: list[Tuple[list[Table] | None, ExitStates]]
    if (n_cores := config.n_cores) is None or n_cores == 1:
        results = execute_one_core_run_ast(config, loading_buffer)
    else:
        results = execute_multicore_ast(
            config, loading_buffer, n_cores)

    # Unpack tables
    outs: list[Table] = [
        sub_table
        for tables, _ in results if tables is not None
        for sub_table in tables
    ]

    # figure out if all workers finished successfully.
    worker_states: list[ExitStates] = [state for _, state in results]
    if all(s == ExitStates.EXIT_SUCCESS for s in worker_states):
        exit_state = ExitStates.EXIT_SUCCESS
    elif any(s == ExitStates.EXIT_SUCCESS for s in worker_states):
        exit_state = ExitStates.EXIT_MIXED
    else:
        exit_state = ExitStates.EXIT_FAIL

    # force finish
    loading_buffer[0] = loading_buffer[1]
    loading.join()

    return outs, exit_state


def top_compile_results(
        outs: list[Table | None], f_name: str,
        old_config: StarBugMainConfig) -> None:
    """
    compiles the results.
    :param outs: the result out tables
    :type outs: list[Table | None]
    :param f_name: the first fits file name.
    :type f_name: str
    :param old_config: the old config
    :type old_config: StarBugMainConfig
    :return: None
    :rtype: None
    """
    config = StarBugMainConfig()
    config.unfreeze()
    config.custom_filter = old_config.custom_filter
    config.fits_images = old_config.fits_images
    config.psf_file_override = old_config.psf_file_override
    config.ast_out_tables = outs
    config._ast_load_psf = False
    config.do_artificial_star_test_results = True
    star_bug_base: StarbugBase = StarbugBase(
        f_name, config, ap_file=config.ap_file,
        bkg_file=config.background_file)
    star_bug_base.run_starbug()


def ast_main(
        argv: list[str], share_memory: SharedMemory,
        loading_buffer: np.ndarray) -> ExitStates:

    config: StarBugMainConfig = ast_parse_argv(argv)

    exit_code: ExitStates
    if config.use_ast_one_time_runs():
        if exit_code := ast_one_time_runs(config):
            share_memory.unlink()
            return exit_code

    print(f"{config.fits_images}")

    if not config.fits_images:
        p_error("must include a fits image to work on\n")
        return ExitStates.EXIT_FAIL

    f_name: str = config.fits_images[0]
    outs: list[Table | None]
    outs, exit_code = execute_ast(config, f_name, loading_buffer)

    # verify it worked
    if exit_code != ExitStates.EXIT_SUCCESS:
        return exit_code

    top_compile_results(outs, f_name, config)

    # Wrapped fix to handle rapid multiprocess teardowns safely
    try:
        share_memory.unlink()
    except FileNotFoundError:
        # The memory handle was already unlinked safely by another thread
        pass
    return exit_code


def ast_main_entry() -> ExitStates:
    """Command line entry point"""
    # globals
    c: np.ndarray = np.array([0, 0, 0], dtype=np.int64)
    share_memory: SharedMemory = (
        shared_memory.SharedMemory(create=True, size=c.nbytes))
    loading_buffer: np.ndarray = np.ndarray(
        c.shape, dtype=c.dtype, buffer=share_memory.buf)
    return ast_main(sys.argv, share_memory, loading_buffer)
