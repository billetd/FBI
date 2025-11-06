from FBI.plotting.plot_main import lompe_scan_plot_vectors, lompe_scan_plot_potential, lompe_scan_plot_potential_polar
from FBI.readwrite import fbi_load_hdf5
import gc
import glob
import re
import os
import datetime as dt


if __name__ == '__main__':

    # Where the plots will be saved
    fbi_dir = '/Users/danielbillett/Data/fbi_dev/test_aligning_records/fitacfs/'

    # List of files to iterate over
    # fbi_files = glob.glob(fbi_dir + "FBI_*.hdf5")
    fbi_file = fbi_dir + 'FBI_20250714050000_20250714050100.hdf5'
    timerange = [dt.datetime(2025, 7, 14, 5, 0),
                 dt.datetime(2025, 7, 14, 5, 1)]

    # for fbi_file in fbi_files:
    if fbi_file:

        # Make a directory to hold the images, if one already doesn't exist
        dirname = fbi_dir + re.search('FBI_(.+?)_', fbi_file).group(1) + '/'
        if not os.path.isdir(dirname):
            os.mkdir(dirname)

        # Read in an FBI hdf5 file
        fbi_data = fbi_load_hdf5(fbi_file, timerange=timerange)

        # Iterate over the records in the file and plot
        for record in fbi_data:
        # for rec in range(7032, 7040):

            lompe_scan_plot_vectors(record, dirname)
            lompe_scan_plot_potential_polar(record, dirname)
            # lompe_scan_plot_potential(record, dirname)
            print(dt.datetime(record['scan_year'][0],
                              record['scan_month'][0],
                              record['scan_day'][0],
                              record['scan_hour'][0],
                              record['scan_minute'][0],
                              record['scan_second'][0]).
                  strftime("%Y-%m-%d %H:%M:%S"))

        del fbi_data
        gc.collect()

