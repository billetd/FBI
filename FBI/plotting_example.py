from FBI.plotting.plot_main import lompe_scan_plot_vectors, lompe_scan_plot_potential, lompe_scan_plot_potential_polar
from FBI.readwrite import fbi_load_hdf5
import gc
import glob
import re
import os
import datetime as dt


if __name__ == '__main__':

    # Where the plots will be saved
    fbi_dir = '/Users/danielbillett/Data/lompe/lompe_test/fitacfs/20231217/'

    # List of files to iterate over
    # fbi_files = glob.glob(fbi_dir + "FBI_*.hdf5")
    fbi_file = fbi_dir + 'FBI_20231217001000_20231217010000.hdf5'
    timerange = [dt.datetime(2023, 12, 17, 0, 10),
                 dt.datetime(2023, 12, 17, 1, 0)]

    # for fbi_file in fbi_files:
    if fbi_file:

        # Make a directory to hold the images, if one already doesn't exist
        dirname = fbi_dir + re.search('FBI_(.+?)_', fbi_file).group(1) + '/'
        if not os.path.isdir(dirname):
            os.mkdir(dirname)

        # Read in an FBI hdf5 file
        print('Reading: ' + fbi_file)
        fbi_data = fbi_load_hdf5(fbi_file, timerange=timerange)

        # Iterate over the records in the file and plot
        for record in fbi_data:
        # for rec in range(7032, 7040):

            lompe_scan_plot_vectors(dirname, record)
            lompe_scan_plot_potential_polar(dirname, record)
            lompe_scan_plot_potential(dirname, record)
            print(dt.datetime(record['scan_year'][0],
                              record['scan_month'][0],
                              record['scan_day'][0],
                              record['scan_hour'][0],
                              record['scan_minute'][0],
                              record['scan_second'][0]).
                  strftime("%Y-%m-%d %H:%M:%S"))

        del fbi_data
        gc.collect()

