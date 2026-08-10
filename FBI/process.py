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
from FBI.fitacf import get_scan_times_widebeam, all_data_make_iterable, median_filter_record, fitacf_get_k_vector_circle
from FBI.fitacf import get_scan_times_old
from pydarn.utils.coordinates import gate2geographic_location
from FBI.grid import lompe_grid_canada
os.environ['RAY_DEDUP_LOGS'] = '0'
os.environ['RAY_BACKEND_LOG_LEVEL'] = 'fatal'
import ray
_worker_apex = None

# Per-process cache of beam/gate geographic positions, keyed by
# (radar_id, max_beams, nrang, rsep, frang). See gate_position_table().
_gate_position_tables = {}


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


    # If no "range_times" is given, work it out based on input data
    # May produce silly scans if there is a mix of normal scan and other modes. Be cautious and only use
    # if you know what data is going in.
    if not range_times:
        # Get scan times within timerange, based on whichever radar started earlier
        # Run the old way if using "scanning data", or the new way if using widebeam data
        if all_data[0][1]['scan'] == 0:
            scan_times, range_times, scan_delta = get_scan_times_old(all_data, timerange)
        else:
            range_times, scan_delta = get_scan_times_widebeam(all_data, timerange)

    # Override scan_delta here to integreate more data per scan
    if scandelta_override is not None:
        scan_delta = scandelta_override

    # Initialise an apexpy object, for magnetic transforms
    apex = apexpy.Apex(range_times[0], refh=300)

    # Get the Superdarn grid
    # This is used for plotting purposes later, e.g. shading darker vectors where there is data
    darn_grid_stuff = grid.sdarn_grid(apex)

    # Remove data we don't use, and make an iterable
    print('Shrinking data...')
    all_data_iterable = all_data_make_iterable(all_data, range_times, scan_delta)

    # Retrive a grid encompassing the SuperDARN Canada PolarDARNs
    canada_grid = lompe_grid_canada(apex)
    del apex # No longer needed

    # Create Emodel object. Pass grid and Hall/Pedersen conductance functions
    model = lompe.Emodel(canada_grid, Hall_Pedersen_conductance=None, ew_regularization_limit=(50, 75))
    del canada_grid  # No longer needed

    # Only initialize Ray if it isn't already running.
    if not ray.is_initialized():
        ray.init(num_cpus=cores, include_dashboard=False, object_store_memory=2 * 1024**3)
        # For debugging. Comment out when not in use
        # ray.init(num_cpus=1, include_dashboard=False, object_store_memory=2 * 1024 ** 3, local_mode=True)

    scan_delta_id       = ray.put(scan_delta)
    darn_grid_stuff_id  = ray.put(darn_grid_stuff)
    med_filter_id       = ray.put(med_filter)
    model_id            = ray.put(model)

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

    # apex = apexpy.Apex(scan_time, refh=300)
    # Initialise apex only once per ray worker and hold on to it.
    # This is because the apxex object can't be serialised with ray.put()
    # Do this minimizes the number of apex intialisations
    global _worker_apex
    if _worker_apex is None:
        _worker_apex = apexpy.Apex(scan_time, refh=300)
        print("Apex initialized on this worker!")

    # Get the data in a format that Lompe likes
    sd_data, rids = prepare_lompe_inputs(_worker_apex, all_data, scan_time, scan_delta, med_filter)

    del all_data # No longer needed

    if sd_data is not None:
        # Run lompe
        try:
            scan_lompe = run_lompe_model(sd_data, model)
        except IndexError:
            scan_lompe = None

        # Collect the model data to save
        if scan_lompe is not None:
            lompe_data = lompe_extract(scan_lompe, _worker_apex, scan_time, darn_grid_stuff, rids)

            # Clean up
            del scan_lompe, sd_data, darn_grid_stuff # No longer needed
            print('Scan complete: ' + scan_time.strftime("%Y-%m-%d %H:%M:%S.%f"))
            return lompe_data

        del sd_data, # No longer needed


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

    # Get data position/value arrays for Lompe
    glat, glon, mlat, mlon, le, ln, le_mag, ln_mag, vlos, vlos_err, rid, ve_mag, vn_mag = (
        get_lompe_data_arrs(apex, all_data, scan_time, scan_delta, med_filter=med_filter))

    coords  = np.vstack((glon, glat))
    los     = np.vstack((le, ln))
    los_mag = np.vstack((le_mag, ln_mag))

    # Make the Lompe data object
    try:
        sd_data = lompe.Data(vlos, coordinates=coords, LOS=los, LOS_mag=los_mag,
                             datatype='convection', error=vlos_err, iweight=1.0)
    except AttributeError:
        print('No data in this scan for some reason. Skipping...')
        return None, None

    return sd_data, rid


def gate_position_table(radar_id, max_beams, nrang, rsep, frang):
    """
    Geographic position of every beam/gate of a radar, as a pair of (max_beams, nrang)
    lookup tables. These depend only on the hdw file and the range parameters, so building
    them once per radar replaces one gate2geographic_location() call per record.

    :param radar_id: pydarn.RadarID
    :param max_beams: number of beams to tabulate
    :param nrang: number of range gates to tabulate
    :param rsep: range separation [km]
    :param frang: distance to the first range gate [km]
    :return: (lats, lons), each shaped (max_beams, nrang)
    """

    nrang = int(nrang)
    key = (radar_id, int(max_beams), nrang, int(rsep), int(frang))
    table = _gate_position_tables.get(key)
    if table is None:
        beams, gates = np.meshgrid(np.arange(max_beams), np.arange(nrang), indexing='ij')
        lat, lon = gate2geographic_location(
            stid=radar_id, beam=beams.ravel(), range_gate=gates.ravel(),
            height=300, center=True, rsep=rsep, frang=frang
        )
        table = (np.asarray(lat).reshape(max_beams, nrang),
                 np.asarray(lon).reshape(max_beams, nrang))
        _gate_position_tables[key] = table

    return table


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
    glon, glat   = [], []
    mlons, mlats = [], []
    vlos, vlos_err = [], []
    le, ln         = [], []
    le_mag, ln_mag = [], []
    ve_mag, vn_mag = [], []
    rid = []

    # For median filtering
    weighting_array = np.array([[[1, 1, 1], [1, 2, 1], [1, 1, 1]],
                                [[2, 2, 2], [2, 4, 2], [2, 2, 2]],
                                [[1, 1, 1], [1, 2, 1], [1, 1, 1]]
                                ])

    for file_index in range(len(all_data)):

        # Station ID
        stid = all_data[file_index][0]['stid']

        # Get position of radar in geographic from hdw files in pyDARN, convert to magnetic
        radar_id = pydarn.RadarID(stid)
        radlat = pydarn.SuperDARNRadars.radars[radar_id].hardware_info.geographic.lat
        radlon = pydarn.SuperDARNRadars.radars[radar_id].hardware_info.geographic.lon
        _rmlat, _rmlon = apex.geo2apex(radlat, radlon, 300)
        radmlat, radmlon = float(_rmlat), float(_rmlon)  # ensure plain floats for math.*

        # Get the indexes for the records which are within half of scan_time
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

        # Collect all valid gate positions and vels in temp arrays
        _lats, _lons   = [], []
        _vlos_signed   = []   # signed: needed so batch function can embed direction in le/ln
        _vlos_mag      = []   # magnitude: what Lompe actually receives
        _vlos_err      = []

        for record in times_in_scan:

            # Ranges with data in it, minus ground scatter
            try:
                slist = all_data[file_index][record]['slist']
            except KeyError:
                continue
            if slist is None:
                continue

            gflg = all_data[file_index][record]['gflg']
            beam = all_data[file_index][record]['bmnum']

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

            # Get coordinates of all beams and gates. These only depend on the hardware and
            # the range parameters, so they come from a per-radar lookup table.
            table_lat, table_lon = gate_position_table(radar_id, max_beams,
                                                       all_data[file_index][record]['nrang'], rsep, frang)
            lat, lon = table_lat[beam, slist], table_lon[beam, slist]

            v    = all_data[file_index][record]['v']
            v_e  = all_data[file_index][record]['v_e']

            # Only keep gates that are not ground scatter, have velocity below 2000m/s, and
            # are above range gate 10. This removes most erroneous data and near-range
            # (E-region) echos.
            keep = (gflg == 0) & (np.abs(v) <= 2000) & (slist > 10)

            # Median filtering. One pass covers every gate of the record.
            if med_filter:
                if not keep.any():
                    continue
                medians, passed = median_filter_record(weighting_array, all_data[file_index],
                                                       record, max_beams)
                vel_range = medians[slist]
                # The median filter can fail if unreliable scatter is found
                keep &= passed[slist]
            else:
                vel_range = v

            # A velocity of exactly zero was dropped by the old `if not vel_range` check,
            # in both the filtered and unfiltered case. Kept for consistency.
            keep &= vel_range != 0.0

            if not keep.any():
                continue

            # Store positions and velocity info
            _lats.extend(lat[keep])
            _lons.extend(lon[keep])
            _vlos_signed.extend(vel_range[keep])
            _vlos_mag.extend(np.abs(vel_range[keep]))
            _vlos_err.extend(np.abs(v_e[keep]))

        if not _lats:
            continue

        # Mag conversion
        _lats_arr = np.array(_lats)
        _lons_arr = np.array(_lons)
        _mlats_arr, _mlons_arr = apex.geo2apex(_lats_arr, _lons_arr, 300)
        _vlos_arr = np.array(_vlos_signed)

        # Get kvectors
        le_arr, ln_arr, le_mag_arr, ln_mag_arr, _, _, ve_mag_arr, vn_mag_arr = \
            fitacf_get_k_vector_circle(
                radlat, radlon, radmlat, radmlon,
                _lats_arr, _lons_arr, _mlats_arr, _mlons_arr, _vlos_arr
            )

        # Add to output lists
        n = len(_lats)
        rid.extend([stid] * n)
        glat.extend(_lats)
        glon.extend(_lons)
        mlats.extend(_mlats_arr.tolist())
        mlons.extend(_mlons_arr.tolist())
        vlos.extend(_vlos_mag)
        vlos_err.extend(_vlos_err)
        le.extend(le_arr.tolist())
        ln.extend(ln_arr.tolist())
        le_mag.extend(le_mag_arr.tolist())
        ln_mag.extend(ln_mag_arr.tolist())
        ve_mag.extend(ve_mag_arr.tolist())
        vn_mag.extend(vn_mag_arr.tolist())

    return (np.array(glat), np.array(glon), np.array(mlats), np.array(mlons),
            np.array(le), np.array(ln), np.array(le_mag), np.array(ln_mag),
            np.array(vlos), np.array(vlos_err), np.array(rid),
            np.array(ve_mag), np.array(vn_mag))


def run_lompe_model(sd_data, model):
    """
    Segregated code to run the lompe model
    """

    # Add all the vectors to the model object
    model.add_data(sd_data)

    # Run inversion
    try:
        model.run_inversion(l1=10, l2=0.1, lapack_driver='gelsy')
    except TypeError:
        # I had the run break on inversion randomly once. Not sure why.
        model = None

    return model