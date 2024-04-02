"""
Example code for processing input fitacf data into lompe outputs
"""

import datetime as dt
import glob
import FBI.fitacf as fitacf
import FBI.process as process
import gc


def simple_process():

    # Locations of files to read
    # Make sure the list of files are only the ones you need, not all your SuperDARN data,
    # otherwise it will read in everything
    fitacf_dir = '/Users/danielbillett/Data/lompe/lompe_test/fitacfs/20231217/imaging/'
    fitacf_files = glob.glob(fitacf_dir+'/*.fitacf')

    # Where to save the lompe outputs
    lompe_dir = '/Users/danielbillett/Data/lompe/lompe_test/fitacfs/20231217/'

    # Times to process between
    # Make sure these times are actually in the fitacfs you have
    start_time = dt.datetime(2023, 12, 17, 0, 10)
    end_time = dt.datetime(2023, 12, 17, 1, 0)

    # Read in 5 at a time (change based on your computers core capacity)
    all_data = fitacf.read_fitacfs(fitacf_files, cores=5)

    # Go and do the rest of the processing
    process.process(all_data, [start_time, end_time], lompe_dir, cores=6, scandelta_override=6)


def dailies_between_dates():
    """
    CAUTION: When reading daily files of the 3.5s data, a LOT of RAM will be consumed
    simply holding the data in memory. Make sure the computer you run this on has enough RAM, else
    swap will be used and processing will be extremely slow.
    Most of the time, this can only be ram on servers with 100+ GB of ram.
    """

    # Locations of the radar directories containing full-day fitacf files,
    # and where to store the FBI output files
    daily_dir = '/sddata/special_experiments/dailies/'
    radar_folders = glob.glob(daily_dir + '*')
    fbi_dir = '/home/billett/data/202401-02_widebeam/fbi_outputs/'

    # Dates to process between
    # We will make one FBI output file for each day
    start_date = dt.datetime(2024, 1, 5)
    end_date = dt.datetime(2024, 2, 14)
    dates = [start_date + i * dt.timedelta(days=1) for i in range((end_date - start_date).days + 1)]

    # Iterate over dates
    for date in dates:
        datestring = date.strftime("%Y%m%d")

        # Check the FBI output files doesn't already exist
        fbi_file_exist = glob.glob(fbi_dir + 'FBI_' + datestring + '*.hdf5')
        if fbi_file_exist:
            print('FBI file for' + datestring + ' already exists. Skipping...')
            continue

        # Find the daily files to read in for each radar
        files = [glob.glob(radar_folder + '/' + datestring + '*.fitacf') for radar_folder in radar_folders]

        # Keep only radars with data
        fitacf_files = [sublist[0] for sublist in files if any(elem is not None for elem in sublist)]

        if not fitacf_files:
            print('No data available for :' + datestring + '. Skipping...')
            continue

        # Read in files
        # CAUTION: This will use a lot of ram for daily 3.5s files
        all_data = fitacf.read_fitacfs(fitacf_files, cores=5)

        # Process this day
        process.process(all_data, [date, date + dt.timedelta(days=1)], fbi_dir,
                        cores=20, scandelta_override=6)

        del all_data
        gc.collect()


simple_process()
# dailies_between_dates()