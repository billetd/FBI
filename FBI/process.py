"""
This module contains code for creating lompe fits between a given timerange, using the already read in
SuperDARN data
"""
import apexpy
import lompe
import datetime as dt
import FBI.grid as grid
import os
import gc
import numpy as np
import warnings
from FBI.readwrite import lompe_extract, fbi_save_hdf5
from FBI.utils import find_indexes_within_time_range, get_kp_iterable
from FBI.fitacf import get_scan_times, all_data_make_iterable, median_filter, fitacf_get_k_vector_circle
from pydarn.utils.coordinates import gate2geographic_location
from FBI.grid import lompe_grid_canada
os.environ['RAY_DEDUP_LOGS'] = '0'
import ray


def process(all_data, timerange, lompe_dir, cores=1, scandelta_override=None):
    """

    :param all_data:
    :param timerange:
    :param lompe_dir:
    :param cores:
    :param scandelta_override:
    """

    # Supress the annoying Pandas warnings because append is depreciated
    warnings.simplefilter(action='ignore', category=FutureWarning)

    # Get scan times within timerange, based on whichever radar started earlier
    scan_times, range_times, scan_delta = get_scan_times(all_data, timerange)

    # Override scan_delta here to integreate more data per scan
    if scandelta_override is not None:
        scan_delta = scandelta_override

    # Initialise an apexpy object, for magnetic transforms
    apex = apexpy.Apex(range_times[0], refh=300)

    # Get Kp indices for time in question, as an interable
    print('Getting Kp indexes...')
    kps = get_kp_iterable(range_times)

    # Get the Superdarn grid
    # This is used for plotting purposes later, e.g. shading darker vectors where there is data
    darn_grid_stuff = grid.sdarn_grid(apex)

    # Remove data we don't use, and make an iterable
    print('Shrinking data...')
    all_data_iterable = all_data_make_iterable(all_data, range_times, scan_delta)

    # lompe_data = lompe_parallel(range_times[0], all_data_iterable[0], kps[0], scan_delta, darn_grid_stuff)

    # Initialise workers for parallelisation (or not) and put in constants
    ray.init(num_cpus=cores)
    scan_delta_id = ray.put(scan_delta)
    darn_grid_stuff_id = ray.put(darn_grid_stuff)
    result_ids = [lompe_parallel.remote(scan_time, this_scan_data, kp, scan_delta_id, darn_grid_stuff_id)
                  for scan_time, this_scan_data, kp
                  in zip(range_times, all_data_iterable, kps)]
    lompes = ray.get(result_ids)
    ray.shutdown()

    # Save as a HDF5 file
    fbi_save_hdf5(lompes, timerange, lompe_dir)


@ray.remote
def lompe_parallel(scan_time, all_data, kp, scan_delta, darn_grid_stuff):
    """
    Code to create a lompe fit for a given scan time. Designed to be paraellelised with ray.
    :param scan_time:
    :param all_data:
    :param kp:
    :param scan_delta:
    :param darn_grid_stuff:
    :return:
    """

    apex = apexpy.Apex(scan_time, refh=300)

    # Get the data in a format that Lompe likes
    sd_data = prepare_lompe_inputs(apex, all_data, scan_time, scan_delta)

    if sd_data is not None:  # Make sure there's data before continuing
        # Run lompe
        scan_lompe = run_lompe_model(scan_time, sd_data, kp)

        # Collect the model data to save
        if scan_lompe is not None:  # I had the run break on inversion randomly once. Not sure why.
            lompe_data = lompe_extract(scan_lompe, apex, scan_time, darn_grid_stuff)

            # Clean up
            del scan_lompe, sd_data, apex
            gc.collect()
            print('Scan complete: ' + scan_time.strftime("%Y-%m-%d %H:%M:%S.%f"))

            return lompe_data


def prepare_lompe_inputs(apex, all_data, scan_time, scan_delta):
    """

    :param apex:
    :param all_data:
    :param scan_time:
    :param scan_delta:
    :return:
    """

    # Get data position/value arrays for Lompe
    glat, glon, mlat, mlon, le, ln, le_mag, ln_mag, vlos, vlos_err, rid, ve_mag, vn_mag = (
        get_lompe_data_arrs(apex, all_data, scan_time, scan_delta, med_filter=True))

    # Make the Lompe data object
    coords = np.vstack((glon, glat))
    los = np.vstack((le, ln))
    los_mag = np.vstack((le_mag, ln_mag))

    try:
        sd_data = lompe.Data(vlos, coordinates=coords, LOS=los, LOS_mag=los_mag, datatype='convection', error=vlos_err,
                             iweight=1.0)
    except AttributeError:
        print('No data in this scan for some reason. Skipping...')
        return None

    return sd_data


def get_lompe_data_arrs(apex, all_data, scan_time, scan_delta, med_filter=False):
    """

    :param apex:
    :param all_data:
    :param scan_time:
    :param scan_delta:
    :param med_filter:
    :return:
    """

    # Arrays that will hold the important parameters
    glon = []
    glat = []
    mlons = []
    mlats = []
    vlos = []
    vlos_err = []
    le = []
    ln = []
    le_mag = []
    ln_mag = []
    ve_geo = []
    vn_geo = []
    ve_mag = []
    vn_mag = []
    rid = []

    for file_index in range(len(all_data)):

        # Station ID
        stid = all_data[file_index][0]['stid']

        # Get the indexes for the records which are within half of of scan_time
        record_times = [dt.datetime(all_data[file_index][x]['time.yr'], all_data[file_index][x]['time.mo'],
                                    all_data[file_index][x]['time.dy'], all_data[file_index][x]['time.hr'],
                                    all_data[file_index][x]['time.mt'], all_data[file_index][x]['time.sc'],
                                    all_data[file_index][x]['time.us'])
                        for x in range(0, len(all_data[file_index]))]
        times_in_scan = find_indexes_within_time_range(record_times, scan_time, catchtime=scan_delta/2)
        max_beams = max([entry["bmnum"] for entry in all_data[file_index]]) + 1
        for record in times_in_scan:

            # Ranges with data in it, minus ground scatter
            try:
                slist = all_data[file_index][record]['slist']
            except KeyError:
                continue
            if slist is None:
                continue
            gflg = all_data[file_index][record]['gflg']

            # Range seperation and frang
            try:
                rsep = all_data[file_index][record]['rsep']
            except KeyError:
                rsep = 45

            # Distance to first range gate
            try:
                frang = all_data[file_index][record]['frang']
            except KeyError:
                frang = 180

            # The current beam
            beam = all_data[file_index][record]['bmnum']

            # Iterate over the gates
            for j, gate in enumerate(slist):
                # Only continue if not ground scatter, velocity is below 2000m/s, and range gates above 10
                # This removes most erroneous data and near-range (E-region) echos
                if gflg[j] == 0 and abs(all_data[file_index][record]['v'][j]) <= 2000 and gate > 10:

                    # Median filtering
                    if med_filter is True:
                        vel_range = median_filter(all_data[file_index], record, max_beams, gate)
                    else:
                        vel_range = all_data[file_index][record]['v'][j]

                    # The median filter can fail if unreliable scatter is found
                    # If so, skip this iteration
                    if not vel_range:
                        continue

                    # Get coordinates of this beam/gate
                    lat, lon = gate2geographic_location(stid=stid, beam=beam, range_gate=gate, height=300,
                                                        center=True, rsep=rsep, frang=frang)
                    mlat, mlon = apex.geo2apex(lat, lon, 300)


                    # Kvectors aren't in fitACF files, so we need to calculate it ourselves
                    # azm = fitacf_get_k_vector(stid, lat, lon, all_data[file_index][record]['v'][j])
                    # Get the unit vectors in east and west directions
                    # le_current, ln_current = los_azimuth2en(azm)
                    (le_current, ln_current, le_mag_current, ln_mag_current, ve_geo_current, vn_geo_current,
                     ve_mag_current, vn_mag_current) = fitacf_get_k_vector_circle(apex, stid, lat, lon, mlat, mlon,
                                                                                  vel_range)

                    # Append to placeholder arrays
                    rid = np.append(rid, all_data[file_index][record]['stid'])
                    glat = np.append(glat, lat)
                    glon = np.append(glon, lon)
                    mlats = np.append(mlats, mlat)
                    mlons = np.append(mlons, mlon)

                    # Needs to be magnitude of the velocity, sign is handled by azimuth
                    vlos = np.append(vlos, abs(vel_range))
                    vlos_err = np.append(vlos_err, abs(all_data[file_index][record]['v_e'][j]))
                    le = np.append(le, le_current)
                    ln = np.append(ln, ln_current)
                    le_mag = np.append(le_mag, le_mag_current)
                    ln_mag = np.append(ln_mag, ln_mag_current)
                    ve_geo = np.append(ve_geo, ve_geo_current)
                    vn_geo = np.append(vn_geo, vn_geo_current)
                    ve_mag = np.append(ve_mag, ve_mag_current)
                    vn_mag = np.append(vn_mag, vn_mag_current)

    return glat, glon, mlats, mlons, le, ln, le_mag, ln_mag, vlos, vlos_err, rid, ve_mag, vn_mag


def run_lompe_model(time, sd_data, kp):
    """

    Segregated code to run the lompe model

    """

    # Retrive a grid encompassing the SuperDARN Canada PolarDARNs
    apex = apexpy.Apex(time, refh=300)
    canada_grid = lompe_grid_canada(apex)

    # Now do conductance using the hardy model
    # Kp = 4  # Manual override
    SH = lompe.conductance.hardy_EUV(canada_grid.lon, canada_grid.lat, int(kp), time, 'hall')
    SP = lompe.conductance.hardy_EUV(canada_grid.lon, canada_grid.lat, int(kp), time, 'pedersen')

    # Create Emodel object. Pass grid and Hall/Pedersen conductance functions
    model = lompe.Emodel(canada_grid, Hall_Pedersen_conductance=(SH, SP))

    # Add all the vectors to the model object
    model.add_data(sd_data)

    # Run inversion
    try:
        model.run_inversion(l1=10, l2=0.1, lapack_driver='gelsy')
        # gtg, ltl = model.run_inversion(l1=10, l2=0)
    except TypeError:
        # I had the run break on inversion randomly once. Not sure why.
        model = None
    del canada_grid

    return model
