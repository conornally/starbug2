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
import copy
import os
from multiprocessing import Pool
from multiprocessing.pool import Pool as PoolType
from typing import Tuple

import numpy as np
from astropy.table import Table

from starbug2.constants import FITS_EXTENSION, ExitStates
from starbug2.core.star_bug_config import StarBugMainConfig
from starbug2.core.starbug_main import StarbugBase
from starbug2.utilities.utils import p_error, split_file_name, printf


def execute_star_bug_main(
        f_name: str, config: StarBugMainConfig) -> StarbugBase | None:
    """
    Worker function to initialise and run standard photometry processes on a
    single file.

    :param f_name: the first fits file name
    :type f_name: str
    :param config: the main config object
    :type config: StarBugMainConfig
    :return: The verified StarbugBase pipeline wrapper instance, or None
             if validation fails
    :rtype: starbug2.StarbugBase or None
    """
    # I've put this here because it takes some time
    from starbug2.core.starbug_main import StarbugBase

    # check file exists
    if not os.path.exists(f_name):
        p_error("can't access %s\n" % f_name)
        return None

    folder, file_name, ext = split_file_name(f_name)

    # check correct extension
    if ext != FITS_EXTENSION:
        p_error("file must be type '.fits' not %s\n" % ext)
        return None

    # extract output files
    ap_file: str | None = config.ap_file
    background_file: str | None = config.background_file

    # find file.
    if config.find_file:
        ap: str = "%s/%s-ap.fits" % (folder, file_name)
        bgd: str = "%s/%s-bgd.fits" % (folder, file_name)
        if os.path.exists(ap) and config.ap_file is None:
            ap_file = ap
        if os.path.exists(bgd) and config.background_file is None:
            background_file = bgd

    # Sorting out the stdout
    if config.verbose_logs:
        printf("-> showing starbug stdout for \"%s\"\n" % f_name)
    elif config.n_cores > 1:
        printf("-> hiding starbug stdout for \"%s\"\n" % f_name)
    else:
        printf("-> %s\n" % f_name)

    # execute
    star_bug_base: StarbugBase | None = StarbugBase(
        f_name, config=config, ap_file=ap_file,
        bkg_file=background_file)
    assert star_bug_base is not None

    result: ExitStates = star_bug_base.run_starbug(config)
    if result != ExitStates.EXIT_SUCCESS:
        return None
    else:
        return star_bug_base


def execute_artificial_stars(
        f_name: str,
        config: StarBugMainConfig) -> Tuple[list[Table] | None, ExitStates]:
    """
    Multiprocessing worker function to run artificial star tests on a given
    file.
    :param f_name: the file to process
    :type f_name: str
    :param config: the config object
    :type config: StarBugMainConfig
    :return: The generated artificial stars recovery catalogue table, or
             None if the file doesn't exist.
    :rtype: Tuple[list[astropy.table.Table] or None., ExitState]
    """
    out: list[Table] | None = None
    exit_state: ExitStates = ExitStates.EXIT_FAIL
    if os.path.exists(f_name):
        star_bug_base: StarbugBase = StarbugBase(
            f_name, config, ap_file=config.ap_file,
            bkg_file=config.background_file, psf=config.ast_psf,
            filter_string=config.custom_filter)
        exit_state = star_bug_base.run_starbug(config)
        out = star_bug_base.ast_test_results
    return out, exit_state


def execute_multicore_ast(
        config: StarBugMainConfig, loading_buffer: np.ndarray,
        n_cores: int) -> list[Tuple[list[Table] | None, ExitStates]]:
    """
    executes tests utilising multicore functionality.
    :param config: the main starbug config
    :type config: StarBugMainConfig
    :param n_cores: the number of cores to utilise
    :type n_cores: int
    :param loading_buffer: the loading buffer
    :type loading_buffer: np.ndarray
    :return: the result tables.
    :rtype: Tuple[list[Table | None], ExitStates]
    """
    # create the initial params
    worker_tasks = [
        (file_name, copy.deepcopy(config))
        for index, file_name in enumerate(config.fits_images)
    ]

    # adjust config to match the instance
    n_cores: int = int(min(n_cores, config.artificial_star_tests_count))
    per_process_n_test: int = int(
        np.ceil(config.artificial_star_tests_count / n_cores))
    per_process_tests_per_save: int = int(
        np.ceil(config.ast_auto_save / n_cores))

    for index, file_name in enumerate(config.fits_images):
        per_process_config: StarBugMainConfig = worker_tasks[index][1]
        per_process_config.unfreeze()
        per_process_config.artificial_star_tests_count = per_process_n_test
        per_process_config.ast_auto_save = per_process_tests_per_save
        per_process_config.verbose_logs = index == 0
        per_process_config.ast_test_index = index
        per_process_config.do_artificial_star_test = True
        per_process_config.ast_loader = loading_buffer
        ast_seed: int | None = config.ast_seed
        if ast_seed is not None:
            per_process_config.ast_seed = (
                ast_seed + (index * per_process_n_test))
        per_process_config.freeze()

    # execute
    pool: PoolType = Pool(processes=n_cores)
    outs = pool.starmap(execute_artificial_stars, worker_tasks)
    pool.close()
    pool.join()

    # return results.
    return outs


def execute_multi_core_main(
        config: StarBugMainConfig) -> list[StarbugBase | None]:
    pool: PoolType = Pool(processes=config.n_cores)

    # this ensures only the first worker executes verbose.
    config.unfreeze()
    config.verbose_logs = False
    config.freeze()

    # create params
    worker_tasks = [
        (file_name, copy.deepcopy(config))
        for index, file_name in enumerate(config.fits_images)
    ]

    # lock the first one to be the only one with verbose logs.
    worker_tasks[0][1].unfreeze()
    worker_tasks[0][1].verbose_logs = True
    worker_tasks[0][1].freeze()

    # trigger the runs.
    starbugs = pool.starmap(execute_star_bug_main, worker_tasks)
    pool.close()
    pool.join()
    return starbugs


def execute_one_core_run_ast(
        config: StarBugMainConfig,
        loading_buffer: np.ndarray) -> (
            list[Tuple[list[Table] | None, ExitStates]]):
    """
    executes tests utilising 1 core.
    :param config: the main starbug config
    :type config: StarBugMainConfig
    :param loading_buffer: the loading buffer
    :type loading_buffer: np.ndarray
    :return: the result tables.
    :rtype: list[Table | None].
    """
    config.unfreeze()
    config.n_cores = 1
    config.do_artificial_star_test = True
    config.ast_loader = loading_buffer
    config.freeze()

    outs = ([execute_artificial_stars(f_name, config)
             for index, f_name in enumerate(config.fits_images)])
    return outs


def execute_one_core_run_main(
        config: StarBugMainConfig) -> list[StarbugBase | None]:
    """
    executes tests utilising 1 core.
    :param config: the main starbug config
    :type config: StarBugMainConfig
    :return: the result tables.
    :rtype: list[StarbugBase | None].
    """
    config.unfreeze()
    config.n_cores = 1
    config.verbose_logs = True
    config.freeze()
    outs: list[StarbugBase | None] = (
        [execute_star_bug_main(f_name, config)
            for index, f_name in enumerate(config.fits_images)])
    return outs
