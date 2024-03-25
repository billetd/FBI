from FBI.plotting.plot_main import lompe_scan_plot_vectors
from FBI.readwrite import fbi_load_hdf5


if __name__ == '__main__':

    # Where the plots will be saved
    save_dir = '/Users/danielbillett/Data/lompe/lompe_test/fitacfs/20231217/'

    # Read in an FBI hdf5 file
    fbi_data = fbi_load_hdf5(save_dir + 'FBI_20231217001000_20231217010000_widebeam.hdf5')

    # Iterate over the records in the file and plot
    for counter, record in enumerate(fbi_data):
        lompe_scan_plot_vectors(save_dir, record)
        print(counter)

