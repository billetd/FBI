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

# Use latex for rendering if installed, fallback if not
_USETEX = bool(shutil.which('latex') and shutil.which('dvipng'))

# font fallback
plt.rcParams['font.sans-serif'] = (['Verdana']
                                   + [font for font in plt.rcParamsDefault['font.sans-serif']
                                      if font != 'Verdana'])

# Image parameters. Increase DPI to increase resolution. 150 = website patch processing. 300 = publication
_DPI = 150
_SUFFIX = '.webp'
_SAVE_KWARGS = {'lossless': True, 'method': 4}

# The start of each kind's image name. extras.plot_fbi_files() relies on these.
PREFIXES = {'vectors': 'vecs_', 'potential': 'pot_', 'potential_polar': 'polar_pot_'}


def scan_time(lompe):
    """
    :param lompe: dict - One record from readwrite.fbi_load_hdf5()
    :return: datetime of the scan
    """

    return dt.datetime(lompe['scan_year'][0], lompe['scan_month'][0], lompe['scan_day'][0], lompe['scan_hour'][0],
                       lompe['scan_minute'][0], lompe['scan_second'][0], lompe['scan_millisec'][0])


def _plot(kind, draw, lompe, path, save, apex):
    """
    The parts common to every kind of plot: skip it if the image already exists, draw it, and save it
    :param kind: str - a key of PREFIXES
    :param draw: function(lompe, time, apex) that draws the plot and returns (fig, ax, ot)
    :param lompe: dict - One record from readwrite.fbi_load_hdf5()
    :param path: str - Directory to save the image in, or None
    :param save: True or False - Save and close the plot, or return it
    :param apex: apexpy.Apex object to use. One is made from the scan time if not given
    :return: (fig, ax, ot, plt) if not saving, otherwise Nones
    """

    time = scan_time(lompe)

    save_path = None
    if path is not None:
        save_path = path + time.strftime(PREFIXES[kind] + "%Y-%m-%d_%H%M%S") + _SUFFIX
        if pathy.isfile(save_path):  # Check plot doesn't already exist
            print('Allready processed: ' + save_path)
            return None

    if apex is None:
        apex = apexpy.Apex(time, refh=300)

    fig, ax, ot = draw(lompe, time, apex)

    if save is True:
        plt.savefig(save_path, dpi=_DPI, bbox_inches='tight', pil_kwargs=_SAVE_KWARGS)
        plt.close('all')
        return None, None, None, None
    else:
        return fig, ax, ot, plt


def _draw_vectors(lompe, time, apex):

    plt.rcParams['text.usetex'] = _USETEX

    # Local axis over Canada
    ax, ot, coord, fig = get_local_axis(apex)
    plt.title(time.strftime("%Y-%m-%d %H:%M:%S"))

    # Line indicating the MLT noon direction
    plot_noon_line(apex, time, coord=coord)

    # Oplot model velocity vectors at grid locations
    plot_vecs_model_darn_grid(lompe, ax, coord=coord)

    return fig, ax, ot


def _draw_potential(lompe, time, apex):

    plt.rcParams['text.usetex'] = _USETEX

    # Local axis over Canada
    ax, ot, coord, fig = get_local_axis(apex)
    plt.title(time.strftime("%Y-%m-%d %H:%M:%S"))

    # Line indicating the MLT noon direction
    plot_noon_line(apex, time, coord=coord)

    # Electric potentials
    plot_potential_contours(lompe, ot, apex, time, coord='mag')

    # Locations of SuperDARN data
    plot_data_locs(lompe, ax, apex=None, time=None, coord=coord)

    return fig, ax, ot


def _draw_potential_polar(lompe, time, apex):

    # Without latex. polplot turns it on when imported, so this has to be set.
    plt.rcParams['text.usetex'] = False

    # Global polar axis
    ax, coord, fig = get_polar_axis(time, apex)
    plt.title(time.strftime("%Y-%m-%d %H:%M:%S"))

    # Electric potentials
    plot_potential_contours(lompe, ax, apex, time, coord='mlt')

    # Locations of SuperDARN data
    plot_data_locs(lompe, ax, apex=apex, time=time, coord=coord)

    # Boundary box
    plot_boundary_box(lompe, ax, apex, time)

    return fig, ax, None


def lompe_scan_plot_vectors(lompe, path=None, save=True, apex=None):
    """
    Fit velocity vectors on the SuperDARN equal area grid, over Canada
    :param lompe: dict - One record from readwrite.fbi_load_hdf5()
    :param path: str - Directory to save the image in
    :param save: True or False - Save and close the plot, or return (fig, ax, ot, plt)
    :param apex: apexpy.Apex object to use. One is made from the scan time if not given
    """

    return _plot('vectors', _draw_vectors, lompe, path, save, apex)


def lompe_scan_plot_potential(lompe, path, save=True, apex=None):
    """
    Electric potential contours and data locations, over Canada
    :param lompe: dict - One record from readwrite.fbi_load_hdf5()
    :param path: str - Directory to save the image in
    :param save: True or False - Save and close the plot, or return (fig, ax, ot, plt)
    :param apex: apexpy.Apex object to use. One is made from the scan time if not given
    """

    return _plot('potential', _draw_potential, lompe, path, save, apex)


def lompe_scan_plot_potential_polar(lompe, path, save=True, apex=None):
    """
    Electric potential contours, data locations and the fit boundary, on a polar plot in MLT
    :param lompe: dict - One record from readwrite.fbi_load_hdf5()
    :param path: str - Directory to save the image in
    :param save: True or False - Save and close the plot, or return (fig, ax, None, plt)
    :param apex: apexpy.Apex object to use. One is made from the scan time if not given
    """

    return _plot('potential_polar', _draw_potential_polar, lompe, path, save, apex)


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
    first = scan_time(records[0]).replace(microsecond=0)
    apex = apexpy.Apex(first, refh=300)

    _shared.update(records=records, path=path, kind=kind, apex=apex)

    # Project the coastlines here, so the workers inherit them
    if kind == 'potential_polar':
        plt.close(get_polar_axis(first, apex)[2])
    else:
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
