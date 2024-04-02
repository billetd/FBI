import apexpy
import matplotlib.pyplot as plt
import datetime as dt
from os import path as pathy
from FBI.plotting.axis import get_local_axis
from FBI.plotting.plot import plot_noon_line, plot_vecs_model_darn_grid, plot_potential_contours, plot_data_locs


def lompe_scan_plot_vectors(path, lompe):
    """

    :param path:
    :param lompe:
    :return:
    """

    plt.rcParams['text.usetex'] = True

    scan_time = dt.datetime(lompe['scan_year'][0], lompe['scan_month'][0], lompe['scan_day'][0], lompe['scan_hour'][0],
                            lompe['scan_minute'][0], lompe['scan_second'][0], lompe['scan_millisec'][0])
    save_path = path + scan_time.strftime("vecs_%Y-%m-%d %H%M%S") + '.png'

    # Check plot doesn't already exist
    if pathy.isfile(save_path) is False:
        # Apex coordinate stuff
        apex = apexpy.Apex(scan_time, refh=300)

        # Local axis over Canada
        ax, ot, coord, fig = get_local_axis(apex)
        plt.title(scan_time.strftime("%Y-%m-%d %H:%M:%S"))

        # Line indicating the MLT noon direction
        plot_noon_line(apex, scan_time, coord=coord)

        # Oplot model velocity vectors at grid locations
        plot_vecs_model_darn_grid(lompe, ax, coord=coord)

        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close('all')
    else:
        print('Allready processed: ' + save_path)


def lompe_scan_plot_potential(path, lompe):
    """

    :param path:
    :param lompe:
    :return:
    """

    plt.rcParams['text.usetex'] = True

    scan_time = dt.datetime(lompe['scan_year'][0], lompe['scan_month'][0], lompe['scan_day'][0], lompe['scan_hour'][0],
                            lompe['scan_minute'][0], lompe['scan_second'][0], lompe['scan_millisec'][0])
    save_path = path + scan_time.strftime("pot_%Y-%m-%d %H%M%S") + '.png'

    # Check plot doesn't already exist
    if pathy.isfile(save_path) is False:
        # Apex coordinate stuff
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

        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close('all')
    else:
        print('Allready processed: ' + save_path)
