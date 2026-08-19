from FBI.plotting.plot_main import plot_records
from FBI.readwrite import fbi_load_hdf5
from FBI.extras import plot_fbi_files
import gc
import glob
import re
import os
import datetime as dt


def plot_everything():
    """
    Turn a whole directory of FBI hdf5 files into plots, skipping whatever is already done.
    Safe to run daily against a directory that keeps growing.
    """

    root = '/Volumes/The Box/FBI/test/batch_test/'

    plot_fbi_files(root + 'fbi_files/', root + 'fbi_plots/', cores=5,
                   kinds=('vectors', 'potential_polar'))


if __name__ == '__main__':

    # # Where the plots will be saved
    # fbi_dir = '/Volumes/The Box/FBI/waves_everywhere/20240213/new_fbi/short/'
    # # fbi_dir = '/Users/danielbillett/Data/FBI/test_data/fitacfs/2025/02/'
    #
    # # List of files to iterate over
    # # fbi_files = glob.glob(fbi_dir + "FBI_*.hdf5")
    # fbi_file = fbi_dir + 'FBI_20240213140000_20240213143000.hdf5'
    # # fbi_file = fbi_dir + 'new_FBI_20250224180000_20250224200000.hdf5'
    # timerange = [dt.datetime(2024, 2, 13, 14, 0),
    #              dt.datetime(2024, 2, 13, 14, 5)]
    #
    # # for fbi_file in fbi_files:
    # if fbi_file:
    #
    #     # Make a directory to hold the images, if one already doesn't exist
    #     dirname = fbi_dir + re.search('FBI_(.+?)_', fbi_file).group(1) + '/'
    #     if not os.path.isdir(dirname):
    #         os.mkdir(dirname)
    #
    #     # Read in an FBI hdf5 file
    #     fbi_data = fbi_load_hdf5(fbi_file, timerange=timerange, as_arrays=True)
    #
    #     # Plot every record. kind can be 'vectors', 'potential' or 'potential_polar'
    #     plot_records(fbi_data, dirname, cores=None, kind='potential_polar')
    #
    #     del fbi_data
    #     gc.collect()

    # Or do a whole directory of FBI files at once, skipping the periods already plotted
    plot_everything()

