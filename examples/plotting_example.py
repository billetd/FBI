"""
Example code for plotting FBI hdf5 files
"""

import os
import re
import datetime as dt
from FBI.plotting.plot_main import plot_records
from FBI.readwrite import fbi_load_hdf5
from FBI.extras import plot_fbi_files


def plot_one_file():
    """
    Plot part of a single FBI hdf5 file, into a directory named after its start time
    """

    fbi_dir = '/Volumes/The Box/FBI/waves_everywhere/20240213/new_fbi/short/'
    fbi_file = fbi_dir + 'FBI_20240213140000_20240213143000.hdf5'
    timerange = [dt.datetime(2024, 2, 13, 14, 0),
                 dt.datetime(2024, 2, 13, 14, 5)]

    # Make a directory to hold the images, if one already doesn't exist
    dirname = fbi_dir + re.search('FBI_(.+?)_', fbi_file).group(1) + '/'
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

    # Read in an FBI hdf5 file
    fbi_data = fbi_load_hdf5(fbi_file, timerange=timerange, as_arrays=True)

    # Plot every record. kind can be 'vectors', 'potential' or 'potential_polar'
    plot_records(fbi_data, dirname, cores=None, kind='potential_polar')


def plot_everything():
    """
    Turn a whole directory of FBI hdf5 files into plots, skipping whatever is already done.
    Safe to run daily against a directory that keeps growing.
    """

    root = '/Volumes/The Box/FBI/test/batch_test/'

    plot_fbi_files(root + 'fbi_files/', root + 'fbi_plots/', cores=5,
                   kinds=('vectors', 'potential_polar'))


if __name__ == '__main__':
    # plot_one_file()
    plot_everything()
