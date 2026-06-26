"""
This module contains code for creating lompe fits between a given timerange, using the already read in
SuperDARN data
"""
import warnings
# Suppress RuntimeWarnings due to not calculating conductances (np.linalg)
warnings.filterwarnings('ignore', category=RuntimeWarning, module='numpy')
import apexpy
import lompe
import datetime as dt
import pydarn
import FBI.grid as grid
import os
import gc
import numpy as np
from FBI.readwrite import lompe_extract, fbi_save_hdf5
from FBI.utils import find_indexes_within_time_range
from FBI.fitacf import get_scan_times_widebeam, all_data_make_iterable, median_filter, fitacf_get_k_vector_circle
from FBI.fitacf import get_scan_times_old
from pydarn.utils.coordinates import gate2geographic_location
from FBI.grid import lompe_grid_canada
os.environ['RAY_DEDUP_LOGS'] = '0'
os.environ['RAY_BACKEND_LOG_LEVEL'] = 'fatal'
import ray


def process(all_data, timerange, lompe_dir, cores=1, med_filter=True, scandelta_override=None, range_times=None):
    """
    :param all_data: list[dict] - List of dictionaries containing fitacf data read in with fitacf.read_fitacfs()
    :param timerange: list[datetime] - Start and end times
    :param lompe_dir: str - Directory to save FBI output file
    :param cores: int - Number of cores to use when multiprocessing. Choose 1 for single core.
    :param med_filter: True or False - Median filter the data before putting into Lompe
    :param scandelta_override: int - Time in seconds to gather data around scans
    :param range_times: list[datetime] - Custom "scan" intervals
    """

    if not range_times:
        if all_data[0][1]['scan'] == 0:
            scan_times, range_times, scan_delta = get_scan_times_old(all_data, timerange)
        else:
            range_times, scan_delta = get_scan_times_widebeam(all_data, timerange)

    if scandelta_override is not None:
        scan_delta = scandelta_override

    # Initialise an apexpy object, for magnetic transforms
    apex = apexpy.Apex(range_times[0], refh=300)
    darn_grid_stuff = grid.sdarn_grid(apex)

    print('Shrinking data...')
    all_data_iterable = all_data_make_iterable(all_data, range_times, scan_delta)

    canada_grid = lompe_grid_canada(apex)
    del apex   # No longer needed

    model = lompe.Emodel(canada_grid, Hall_Pedersen_conductance=None, ew_regularization_limit=(50, 75))
    del canada_grid  # No longer needed

    # Only initialize Ray if it isn't already running.
    if not ray.is_initialized():
        ray.init(num_cpus=cores, include_dashboard=False)

    scan_delta_id     = ray.put(scan_delta)
    darn_grid_stuff_id = ray.put(darn_grid_stuff)
    med_filter_id     = ray.put(med_filter)
    model_id          = ray.put(model)

    # Bounded task submission via ray.wait().
    # keep at most `cores` tasks in-flight at any time, submitting the
    # next one only when a slot frees up.
    all_pairs = list(zip(range_times, all_data_iterable))

    del darn_grid_stuff, model
    del all_data_iterable
    del all_data # No longer needed
    gc.collect()

    n_total = len(all_pairs)
    lompes  = [None] * n_total

    # Each entry: (result_id, original_index)
    pending   = []
    submitted = 0

    # Submit the first batch (up to `cores` tasks)
    for _ in range(min(cores, n_total)):
        scan_time, this_scan_data = all_pairs[submitted]
        rid = lompe_parallel.remote(
            scan_time, this_scan_data,
            scan_delta_id, darn_grid_stuff_id, med_filter_id, model_id
        )
        pending.append((rid, submitted))
        submitted += 1

    # Rolling window: as each task finishes, collect its result and launch the next
    while pending:
        pending_ids = [p[0] for p in pending]
        done_ids, _ = ray.wait(pending_ids, num_returns=1, timeout=600)

        if not done_ids:
            # A task is taking longer than 1 minute — unusual, so warn the user. Maybe too long a scandeltaoverride?
            print("Warning: task is taking unusually long, still waiting. Is your scan_delta too long?")
            continue

        done_id = done_ids[0]
        original_idx = next(idx for rid, idx in pending if rid == done_id)
        pending = [(rid, idx) for rid, idx in pending if rid != done_id]

        lompes[original_idx] = ray.get(done_id)

        # Submit the next pending task now that a worker slot has freed up
        if submitted < n_total:
            scan_time, this_scan_data = all_pairs[submitted]
            rid = lompe_parallel.remote(
                scan_time, this_scan_data,
                scan_delta_id, darn_grid_stuff_id, med_filter_id, model_id
            )
            pending.append((rid, submitted))
            submitted += 1

    ray.shutdown()

    fbi_save_hdf5(lompes, timerange, lompe_dir)

@ray.remote
def lompe_parallel(scan_time, all_data, scan_delta, darn_grid_stuff, med_filter, model):
    """
    Code to create a lompe fit for a given scan time. Designed to be parallelised with ray.
    :param scan_time:
    :param all_data:
    :param scan_delta:
    :param darn_grid_stuff:
    :param med_filter:
    :param model:
    :return:
    """

    apex = apexpy.Apex(scan_time, refh=300)

    # Extract compact numpy arrays from the raw scan records
    sd_data, rids = prepare_lompe_inputs(apex, all_data, scan_time, scan_delta, med_filter)

    del all_data # No longer needed

    if sd_data is not None:
        try:
            scan_lompe = run_lompe_model(sd_data, model)
        except IndexError:
            scan_lompe = None

        if scan_lompe is not None:
            lompe_data = lompe_extract(scan_lompe, apex, scan_time, darn_grid_stuff, rids)

            del scan_lompe, sd_data, apex, darn_grid_stuff # No longer needed
            print('Scan complete: ' + scan_time.strftime("%Y-%m-%d %H:%M:%S.%f"))
            return lompe_data

        del sd_data, apex # No longer needed


    return None


def prepare_lompe_inputs(apex, all_data, scan_time, scan_delta, med_filter):
    """
    :param apex:
    :param all_data:
    :param scan_time:
    :param scan_delta:
    :param med_filter: True or false
    :return:
    """

    glat, glon, mlat, mlon, le, ln, le_mag, ln_mag, vlos, vlos_err, rid, ve_mag, vn_mag = (
        get_lompe_data_arrs(apex, all_data, scan_time, scan_delta, med_filter=med_filter))

    coords  = np.vstack((glon, glat))
    los     = np.vstack((le, ln))
    los_mag = np.vstack((le_mag, ln_mag))

    try:
        sd_data = lompe.Data(vlos, coordinates=coords, LOS=los, LOS_mag=los_mag,
                             datatype='convection', error=vlos_err, iweight=1.0)
    except AttributeError:
        print('No data in this scan for some reason. Skipping...')
        return None, None

    return sd_data, rid


def get_lompe_data_arrs(apex, all_data, scan_time, scan_delta, med_filter=False):
    """
    :param apex:
    :param all_data:
    :param scan_time:
    :param scan_delta:
    :param med_filter:
    :return:
    """

    glon     = []
    glat     = []
    mlons    = []
    mlats    = []
    vlos     = []
    vlos_err = []
    le       = []
    ln       = []
    le_mag   = []
    ln_mag   = []
    ve_mag   = []
    vn_mag   = []
    rid      = []

    for file_index in range(len(all_data)):

        stid = all_data[file_index][0]['stid']

        record_times = [
            dt.datetime(
                all_data[file_index][x]['time.yr'], all_data[file_index][x]['time.mo'],
                all_data[file_index][x]['time.dy'], all_data[file_index][x]['time.hr'],
                all_data[file_index][x]['time.mt'], all_data[file_index][x]['time.sc'],
                all_data[file_index][x]['time.us']
            )
            for x in range(len(all_data[file_index]))
        ]
        times_in_scan = find_indexes_within_time_range(record_times, scan_time, catchtime=scan_delta / 2)
        max_beams = max([entry["bmnum"] for entry in all_data[file_index]]) + 1

        for record in times_in_scan:

            try:
                slist = all_data[file_index][record]['slist']
            except KeyError:
                continue
            if slist is None:
                continue

            gflg = all_data[file_index][record]['gflg']

            try:
                rsep = all_data[file_index][record]['rsep']
            except KeyError:
                rsep = 45

            try:
                frang = all_data[file_index][record]['frang']
            except KeyError:
                frang = 180

            beam = all_data[file_index][record]['bmnum']

            for j, gate in enumerate(slist):
                if gflg[j] == 0 and abs(all_data[file_index][record]['v'][j]) <= 2000 and gate > 10:

                    if med_filter is True:
                        vel_range = median_filter(all_data[file_index], record, max_beams, gate)
                    else:
                        vel_range = all_data[file_index][record]['v'][j]

                    if not vel_range:
                        continue

                    lat, lon = gate2geographic_location(
                        stid=pydarn.RadarID(stid), beam=beam, range_gate=gate,
                        height=300, center=True, rsep=rsep, frang=frang
                    )
                    mlat, mlon = apex.geo2apex(lat, lon, 300)

                    (le_current, ln_current, le_mag_current, ln_mag_current,
                     ve_geo_current, vn_geo_current,  # returned by the function but not stored
                     ve_mag_current, vn_mag_current) = fitacf_get_k_vector_circle(
                        apex, stid, lat, lon, mlat, mlon, vel_range
                    )

                    rid.append(all_data[file_index][record]['stid'])
                    glat.append(lat)
                    glon.append(lon)
                    mlats.append(mlat)
                    mlons.append(mlon)
                    vlos.append(abs(vel_range))
                    vlos_err.append(abs(all_data[file_index][record]['v_e'][j]))
                    le.append(le_current)
                    ln.append(ln_current)
                    le_mag.append(le_mag_current)
                    ln_mag.append(ln_mag_current)
                    ve_mag.append(ve_mag_current)
                    vn_mag.append(vn_mag_current)

    return (np.array(glat), np.array(glon), np.array(mlats), np.array(mlons),
            np.array(le), np.array(ln), np.array(le_mag), np.array(ln_mag),
            np.array(vlos), np.array(vlos_err), np.array(rid),
            np.array(ve_mag), np.array(vn_mag))


def run_lompe_model(sd_data, model):
    """
    Segregated code to run the lompe model
    """
    model.add_data(sd_data)

    try:
        model.run_inversion(l1=10, l2=0.1, lapack_driver='gelsy')
    except TypeError:
        model = None

    return model