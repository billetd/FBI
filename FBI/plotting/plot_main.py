"""
Top-level functions for producing complete FBI scan plots: build an axis, draw the relevant
elements from FBI.plotting.plot onto it, and save (or return) the figure.
See FBI/plotting_example.py for typical usage.
"""
import apexpy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
import os

from FBI.readwrite import record_datetime
from FBI.plotting.axis import get_local_axis, get_polar_axis
from FBI.plotting.plot import plot_noon_line, plot_vecs_model_darn_grid, plot_potential_contours, plot_data_locs, \
    plot_boundary_box

DEFAULT_REFH_KM = 300
DEFAULT_DPI = 300


def _build_save_path(directory, prefix, scan_time):
    if directory is None:
        return None
    filename = scan_time.strftime(f"{prefix}_%Y-%m-%d %H%M%S") + '.png'
    return os.path.join(directory, filename)


def _safe_tight_bbox(fig, pad_inches=0.1):
    """
    A `bbox_inches='tight'` replacement, robust to a cartopy bug where a GeoAxes with
    `gridlines(draw_labels=True)` can report a tight bbox with NaN/inf extents. Left unguarded,
    that NaN poisons matplotlib's usual tight-bbox union and silently crops saved figures down to
    whatever other (finite) artist happens to be on the axes - e.g. just the colourbar. Any axis
    that reports a non-finite tight bbox falls back to its plain window extent instead.
    """

    renderer = fig.canvas.get_renderer()
    bboxes = []
    for ax in fig.get_axes():
        bbox = ax.get_tightbbox(renderer)
        if not np.all(np.isfinite(bbox.extents)):
            # The window extent excludes artists like a title that sit outside the axes frame,
            # so add those back in individually where they're safe to do so.
            bbox = ax.get_window_extent(renderer)
            for text in (*ax.texts, ax.title):
                if text.get_visible() and text.get_text():
                    text_bbox = text.get_window_extent(renderer)
                    if np.all(np.isfinite(text_bbox.extents)):
                        bbox = Bbox.union([bbox, text_bbox])
        bboxes.append(bbox)

    tight = Bbox.union(bboxes).transformed(fig.dpi_scale_trans.inverted())
    return tight.padded(pad_inches)


def _set_map_title(ax, text):
    """
    Place a title above `ax` as a plain text label rather than via ax.set_title(). set_title()
    relies on matplotlib's automatic title-positioning, which - on a GeoAxes with
    gridlines(draw_labels=True) - inherits a NaN bbox from a cartopy Gridliner bug and breaks
    bbox_inches='tight' entirely (see _safe_tight_bbox).
    """
    ax.text(0.5, 1.02, text, transform=ax.transAxes, ha='center', va='bottom')


def _finalize(fig, save, save_path, dpi=DEFAULT_DPI):
    if not save:
        return fig
    fig.savefig(save_path, dpi=dpi, bbox_inches=_safe_tight_bbox(fig))
    plt.close(fig)
    return None


def lompe_scan_plot_vectors(record, path=None, save=True, use_latex=True, refh=DEFAULT_REFH_KM):
    """
    Plot the Lompe fit velocity field over the SuperDARN Canada sector.
    :param record: scan dict, as returned by FBI.readwrite.fbi_load_hdf5().
    :param path: Directory to save the plot into. If a plot for this scan already exists there,
        plotting is skipped and None is returned. If None, nothing is checked/saved on disk (see `save`).
    :param save: If True, save the figure and close it, returning None. If False, leave the figure
        open and return it.
    :param use_latex: Render text with LaTeX (requires a working LaTeX install).
    :param refh: Apex reference height in km, used for the noon line and coastlines.
    :return: The matplotlib Figure if save is False, else None.
    """

    scan_time = record_datetime(record)
    save_path = _build_save_path(path, 'vecs', scan_time)
    if save_path is not None and os.path.isfile(save_path):
        print('Already processed: ' + save_path)
        return None

    plt.rcParams['text.usetex'] = use_latex

    apex = apexpy.Apex(scan_time, refh=refh)
    ax, projection, coord, fig = get_local_axis(apex, refh=refh)
    _set_map_title(ax, scan_time.strftime("%Y-%m-%d %H:%M:%S"))

    plot_noon_line(apex, scan_time, ax, coord=coord)
    plot_vecs_model_darn_grid(record, ax, coord=coord)

    return _finalize(fig, save, save_path)


def lompe_scan_plot_potential(record, path=None, save=True, use_latex=True, refh=DEFAULT_REFH_KM):
    """
    Plot the Lompe fit electric potential over the SuperDARN Canada sector.
    :param record: scan dict, as returned by FBI.readwrite.fbi_load_hdf5().
    :param path: Directory to save the plot into. If a plot for this scan already exists there,
        plotting is skipped and None is returned. If None, nothing is checked/saved on disk (see `save`).
    :param save: If True, save the figure and close it, returning None. If False, leave the figure
        open and return it.
    :param use_latex: Render text with LaTeX (requires a working LaTeX install).
    :param refh: Apex reference height in km, used for the noon line and coastlines.
    :return: The matplotlib Figure if save is False, else None.
    """

    scan_time = record_datetime(record)
    save_path = _build_save_path(path, 'pot', scan_time)
    if save_path is not None and os.path.isfile(save_path):
        print('Already processed: ' + save_path)
        return None

    plt.rcParams['text.usetex'] = use_latex

    apex = apexpy.Apex(scan_time, refh=refh)
    ax, projection, coord, fig = get_local_axis(apex, refh=refh)
    _set_map_title(ax, scan_time.strftime("%Y-%m-%d %H:%M:%S"))

    plot_noon_line(apex, scan_time, ax, coord=coord)
    plot_potential_contours(record, ax, apex, scan_time, coord=coord, projection=projection)
    plot_data_locs(record, ax, coord=coord)

    return _finalize(fig, save, save_path)


def lompe_scan_plot_potential_polar(record, path=None, save=True, use_latex=True, refh=DEFAULT_REFH_KM):
    """
    Plot the Lompe fit electric potential on a full magnetic-local-time polar map.
    :param record: scan dict, as returned by FBI.readwrite.fbi_load_hdf5().
    :param path: Directory to save the plot into. If a plot for this scan already exists there,
        plotting is skipped and None is returned. If None, nothing is checked/saved on disk (see `save`).
    :param save: If True, save the figure and close it, returning None. If False, leave the figure
        open and return it.
    :param use_latex: Render text with LaTeX (requires a working LaTeX install).
    :param refh: Apex reference height in km.
    :return: The matplotlib Figure if save is False, else None.
    """

    scan_time = record_datetime(record)
    save_path = _build_save_path(path, 'polar_pot', scan_time)
    if save_path is not None and os.path.isfile(save_path):
        print('Already processed: ' + save_path)
        return None

    plt.rcParams['text.usetex'] = use_latex

    apex = apexpy.Apex(scan_time, refh=refh)
    pax, coord, fig = get_polar_axis(scan_time, apex)
    pax.ax.set_title(scan_time.strftime("%Y-%m-%d %H:%M:%S"))

    plot_potential_contours(record, pax, apex, scan_time, coord=coord)
    plot_data_locs(record, pax, apex=apex, time=scan_time, coord=coord)
    plot_boundary_box(record, pax, apex, scan_time)

    return _finalize(fig, save, save_path)
