import apexpy
import matplotlib.pyplot as plt
import datetime as dt
import shutil
import time
from os import path as pathy
from FBI.parallel import resolve_cores, forked_pool, bounded_imap, report_progress
from FBI.plotting.axis import get_local_axis, get_polar_axis
from FBI.plotting.plot import plot_noon_line, plot_vecs_model_darn_grid, plot_potential_contours, plot_data_locs, \
    plot_boundary_box

# Use latex for rendering if install, fallback if not
_USETEX = bool(shutil.which('latex') and shutil.which('dvipng'))

# font fallback
plt.rcParams['font.sans-serif'] = (['Verdana']
                                   + [font for font in plt.rcParamsDefault['font.sans-serif']
                                      if font != 'Verdana'])

def lompe_scan_plot_vectors(lompe, path=None, save=True, apex=None):
    """

    :param path:
    :param lompe:
    :param save:
    :param apex: apexpy.Apex object to use. One is made from the scan time if not given
    :return:
    """

    plt.rcParams['text.usetex'] = _USETEX
    scan_time = dt.datetime(lompe['scan_year'][0], lompe['scan_month'][0], lompe['scan_day'][0], lompe['scan_hour'][0],
                            lompe['scan_minute'][0], lompe['scan_second'][0], lompe['scan_millisec'][0])

    if path is not None:
        save_path = path + scan_time.strftime("vecs_%Y-%m-%d_%H%M%S") + '.png'
        if pathy.isfile(save_path) is False:  # Check plot doesn't already exist
            go = True
        else:
            go = False
    else:
        go = True

    if go is True:
        # Apex coordinate stuff
        if apex is None:
            apex = apexpy.Apex(scan_time, refh=300)

        # Local axis over Canada
        ax, ot, coord, fig = get_local_axis(apex)
        plt.title(scan_time.strftime("%Y-%m-%d %H:%M:%S"))

        # Line indicating the MLT noon direction
        plot_noon_line(apex, scan_time, coord=coord)

        # Oplot model velocity vectors at grid locations
        plot_vecs_model_darn_grid(lompe, ax, coord=coord)

        if save is True:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            plt.close('all')
            return None, None, None, None
        else:
            return fig, ax, ot, plt
    else:
        print('Allready processed: ' + save_path)


def lompe_scan_plot_potential(lompe, path, save=True, apex=None):
    """

    :param path:
    :param lompe:
    :param save:
    :param apex: apexpy.Apex object to use. One is made from the scan time if not given
    :return:
    """

    plt.rcParams['text.usetex'] = _USETEX

    scan_time = dt.datetime(lompe['scan_year'][0], lompe['scan_month'][0], lompe['scan_day'][0], lompe['scan_hour'][0],
                            lompe['scan_minute'][0], lompe['scan_second'][0], lompe['scan_millisec'][0])
    if path is not None:
        save_path = path + scan_time.strftime("pot_%Y-%m-%d_%H%M%S") + '.png'
        if pathy.isfile(save_path) is False:  # Check plot doesn't already exist
            go = True
        else:
            go = False
    else:
        go = True

    if go is True:
        # Apex coordinate stuff
        if apex is None:
            apex = apexpy.Apex(scan_time, refh=300)

        # Local axis over Canada
        ax, ot, coord, fig = get_local_axis(apex)
        plt.title(scan_time.strftime("%Y-%m-%d %H:%M:%S"))

        # Line indicating the MLT noon direction
        plot_noon_line(apex, scan_time, coord=coord)

        # Electric potentials
        plot_potential_contours(lompe, ot, apex, scan_time, coord='mag')

        # Locations of SuperDARN data
        plot_data_locs(lompe, ax, apex=None, time=None, coord=coord)

        if save is True:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            plt.close('all')
            return None, None, None, None
        else:
            return fig, ax, ot, plt
    else:
        print('Allready processed: ' + save_path)


def lompe_scan_plot_potential_polar(lompe, path, save=True, apex=None):
    """

    :param path:
    :param lompe:
    :param save:
    :param apex: apexpy.Apex object to use. One is made from the scan time if not given
    :return:
    """

    # plt.rcParams['text.usetex'] = True

    scan_time = dt.datetime(lompe['scan_year'][0], lompe['scan_month'][0], lompe['scan_day'][0], lompe['scan_hour'][0],
                            lompe['scan_minute'][0], lompe['scan_second'][0], lompe['scan_millisec'][0])

    if path is not None:
        save_path = path + scan_time.strftime("polar_pot_%Y-%m-%d_%H%M%S") + '.png'
        if pathy.isfile(save_path) is False:  # Check plot doesn't already exist
            go = True
        else:
            go = False
    else:
        go = True

    if go is True:
        # Apex coordinate stuff
        if apex is None:
            apex = apexpy.Apex(scan_time, refh=300)

        # Global polar axis
        ax, coord, fig = get_polar_axis(scan_time, apex)
        plt.title(scan_time.strftime("%Y-%m-%d %H:%M:%S"))

        # Electric potentials
        plot_potential_contours(lompe, ax, apex, scan_time, coord='mlt')

        # Locations of SuperDARN data
        plot_data_locs(lompe, ax, apex=apex, time=scan_time, coord=coord)

        # Boundary box
        plot_boundary_box(lompe, ax, apex, scan_time)

        if save is True:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            plt.close('all')
            return None, None, None, None
        else:
            return fig, ax, None, plt
    else:
        print('Allready processed: ' + save_path)


_PLOT_FUNCS = {'vectors': lompe_scan_plot_vectors,
               'potential': lompe_scan_plot_potential,
               'potential_polar': lompe_scan_plot_potential_polar}

# What the workers plot from. Filled in before forking, so a task is just a record index.
_shared = {}


def _plot_one(index):
    """
    Plot a single record. This is what each worker runs.
    :param index: int - position of the record in the list given to plot_records()
    """

    _PLOT_FUNCS[_shared['kind']](_shared['records'][index], _shared['path'],
                                 apex=_shared['apex'])
    plt.close('all')


def plot_records(records, path, cores=None, kind='vectors'):
    """
    Plot a set of records, one image each
    :param records: list[dict] - records from readwrite.fbi_load_hdf5()
    :param path: str - Directory to save the images in
    :param cores: int - Number of worker processes. None uses every CPU available. Choose 1
                  to plot in this process with no pool.
    :param kind: str - 'vectors', 'potential' or 'potential_polar'
    """

    if kind not in _PLOT_FUNCS:
        raise ValueError('kind must be one of ' + str(sorted(_PLOT_FUNCS)) + ', not ' + repr(kind))

    n_total = len(records)
    if not n_total:
        return
    cores = min(resolve_cores(cores), n_total)

    # One apex for the whole run. Every record in a file shares an epoch.
    first = records[0]
    apex = apexpy.Apex(dt.datetime(first['scan_year'][0], first['scan_month'][0], first['scan_day'][0],
                                   first['scan_hour'][0], first['scan_minute'][0], first['scan_second'][0]),
                       refh=300)

    _shared.update(records=records, path=path, kind=kind, apex=apex)

    # Project the coastlines here, so the workers inherit them
    if kind != 'potential_polar':
        plt.close(get_local_axis(apex)[3])

    try:
        started = time.monotonic()
        if cores == 1:
            # Single process, so exceptions and profilers behave normally
            for index in range(n_total):
                _plot_one(index)
                report_progress(index + 1, n_total, started, unit='plots')
        else:
            # One record per task. Each worker saves its own image.
            with forked_pool(cores) as pool:
                for done, _ in enumerate(bounded_imap(pool, _plot_one, n_total, 2 * cores), 1):
                    report_progress(done, n_total, started, unit='plots')
        print()
    finally:
        _shared.clear()
