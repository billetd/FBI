from FBI.plotting.plot_main import lompe_scan_plot_vectors
from FBI.readwrite import fbi_load_hdf5
import gc
import glob
import re
import os


if __name__ == '__main__':

    # Where the plots will be saved
    fbi_dir = '/Users/danielbillett/Data/lompe/2024_widebeam_h5s/'

    # List of files to iterate over
    # fbi_files = glob.glob(fbi_dir + "FBI_*.hdf5")
    fbi_file = fbi_dir + 'FBI_20240111000000_20240112000000.hdf5'

    # for fbi_file in fbi_files:
    if fbi_file:

        # Make a directory to hold the images, if one already doesn't exist
        dirname = fbi_dir + re.search('FBI_(.+?)_', fbi_file).group(1) + '/'
        if not os.path.isdir(dirname):
            os.mkdir(dirname)

        # Read in an FBI hdf5 file
        print('Reading: ' + fbi_file)
        fbi_data = fbi_load_hdf5(fbi_file)

        # Iterate over the records in the file and plot
        # for counter, record in enumerate(fbi_data):
        for rec in range(7032, 7035):

            lompe_scan_plot_vectors(dirname, fbi_data[rec])
            print(fbi_file, rec)

        del fbi_data
        gc.collect()

