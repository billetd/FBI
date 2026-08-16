"""
This module contains code for spreading work across the cores of a single machine
"""

import gc
import contextlib
import multiprocessing
import os
import threading
import warnings
from collections import deque


def resolve_cores(cores=None):
    """
    Work out how many worker processes to use
    :param cores: int or None. None decides based on the machine
    :return: int, at least 1
    """

    if cores is not None:
        return max(1, int(cores))

    env = os.environ.get('FBI_CORES')
    if env:
        return max(1, int(env))

    try:
        # Number of CPUs this process is allowed on, not the number in the machine.
        # Servers usually pin a job to a subset with cgroups or Slurm.
        return len(os.sched_getaffinity(0))
    except AttributeError:  # Not on linux
        return os.cpu_count() or 1


@contextlib.contextmanager
def forked_pool(cores):
    """
    A pool of workers which inherit everything this process has already built, rather
    than being sent a copy of it
    :param cores: int - number of worker processes
    :return: context manager yielding a multiprocessing.Pool
    """

    if 'fork' not in multiprocessing.get_all_start_methods():
        raise RuntimeError(
            'FBI needs the fork start method, which this platform does not provide. '
            'Run with cores=1 to process in a single process instead.'
        )

    # Only the calling thread survives a fork, so a lock held by any other thread stays
    # locked forever in the workers. Nothing in FBI makes threads, but a library might.
    if threading.active_count() > 1:
        warnings.warn(
            f'Forking a worker pool while {threading.active_count()} threads are '
            f'running. This is not safe; if the workers hang, that is why.',
            RuntimeWarning, stacklevel=3,
        )

    # A collection in a worker walks every inherited object to update its generation,
    # which copies the pages we are trying to share. Freezing puts everything that
    # already exists somewhere the collector never looks.
    gc.collect()
    gc.freeze()
    pool = multiprocessing.get_context('fork').Pool(processes=cores)
    try:
        yield pool
    finally:
        pool.terminate()
        pool.join()
        gc.unfreeze()


def bounded_imap(pool, func, n_items, window):
    """
    Map func over range(n_items), yielding (index, result) in order, with no more than
    `window` tasks in flight. Pool.imap() buffers every result that comes back, so it
    piles up if the consumer is slower than the workers.
    :param pool: multiprocessing.Pool
    :param func: module level function taking a single int
    :param n_items: int - number of items
    :param window: int - maximum tasks submitted but not yet collected
    :return: generator of (index, result)
    """

    inflight = deque()
    submitted = 0

    while submitted < n_items or inflight:
        while submitted < n_items and len(inflight) < window:
            inflight.append((submitted, pool.apply_async(func, (submitted,))))
            submitted += 1

        index, pending = inflight.popleft()
        yield index, pending.get()
