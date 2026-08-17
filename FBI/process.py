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
import FBI.readwrite as readwrite
import gc
import os
import time
import numpy as np
from FBI.readwrite import lompe_extract, FBIWriter
from FBI.parallel import resolve_cores, forked_pool, bounded_imap, report_progress
from FBI.utils import find_indexes_within_time_range
from FBI.fitacf import get_scan_times_widebeam, all_data_make_iterable, median_filter_record, fitacf_get_k_vector_circle
from FBI.fitacf import get_scan_times_old
from pydarn.utils.coordinates import gate2geographic_location
from FBI.grid import lompe_grid_canada

# apexpy object used for every scan. An Apex object can't be pickled, so it is built here
# and the workers inherit it.
_worker_apex = None

# Per-process cache of beam/gate geographic positions, keyed by
# (radar_id, max_beams, nrang, rsep, frang). See gate_position_table().
_gate_position_tables = {}

# The model, grids and record windows the workers read. Filled in before forking, so the
# workers inherit it and a task is just a scan index.
_shared = {}


def process(all_data, timerange, lompe_dir, cores=None, med_filter=True, scandelta_override=None, range_times=None,
            cache_geometry=True):
    """
    :param all_data: list[dict] - List of dictionaries containing fitacf data read in with fitacf.read_fitacfs()
    :param timerange: list[datetime] - Start and end times
    :param lompe_dir: str - Directory to save FBI output file
    :param cores: int - Number of worker processes. None uses every CPU available. Choose 1
                  to fit in this process with no pool, for profiling or debugging.
    :param med_filter: True or False - Median filter the data before putting into Lompe
    :param scandelta_override: int - Time in seconds to gather data around scans
    :param range_times: list[datetime] - Custom "scan" intervals
    :param cache_geometry: True or False - Reuse the SECS/apex matrices that depend only on the
                           grids between scans. Much faster, and costs a few hundred MB for
                           the run rather than per core.
    """

    cores = resolve_cores(cores)

    # Clear anything left by a previous call, e.g. from extras.process_date()
    _reset_caches()

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

    # Build the data-density grid here rather than in every worker's first inversion
    model.prepare_biggrid()

    # Build everything the workers need before forking, so they share it rather than
    # each being sent a copy
    _prime_shared_state(range_times, all_data_iterable, scan_delta, darn_grid_stuff,
                        med_filter, model, cache_geometry)

    n_total = len(range_times)
    cores = min(cores, max(1, n_total))  # No point in more workers than scans

    del darn_grid_stuff, model
    del all_data_iterable
    del all_data # No longer needed
    gc.collect()

    try:
        started = time.monotonic()
        with FBIWriter(timerange, lompe_dir) as writer:
            if cores == 1:
                # Single process, so exceptions and profilers behave normally
                for index in range(n_total):
                    writer.write(index, _lompe_one_scan(index))
                    report_progress(index + 1, n_total, started)
            else:
                # One scan per task, for the best load balancing. The window limits how
                # many finished scans can be waiting to be written.
                with forked_pool(cores) as pool:
                    for index, result in bounded_imap(pool, _lompe_one_scan, n_total, 2 * cores):
                        writer.write(index, result)
                        report_progress(index + 1, n_total, started)
            print('\nWrote ' + str(writer.n_written) + ' of ' + str(n_total) + ' scans')
    finally:
        # Drop the model and record windows before the caller moves on to the next chunk
        _reset_caches()
        gc.collect()


def _prime_shared_state(range_times, all_data_iterable, scan_delta, darn_grid_stuff,
                        med_filter, model, cache_geometry):
    """
    Set up the state the workers read, before the pool is forked
    :param range_times: list[datetime] - scan times
    :param all_data_iterable: list - the record window for each scan
    :param scan_delta: int - seconds of data to gather around each scan
    :param darn_grid_stuff: dict from FBI.grid.sdarn_grid()
    :param med_filter: True or False
    :param model: lompe Emodel
    :param cache_geometry: True or False
    """

    global _worker_apex

    # Must come after the Emodel is built. Emodel makes its own Apex at epoch 2015, and
    # apexpy holds the epoch in Fortran state shared by every Apex in the process, so
    # building ours last puts the epoch back to the one the data wants.
    _worker_apex = apexpy.Apex(range_times[0], refh=300)

    if cache_geometry:
        readwrite.prime_geometry_cache(model, _worker_apex, darn_grid_stuff)

    _freeze_model_arrays(model)

    _shared.update(range_times=range_times,
                   all_data_iterable=all_data_iterable,
                   scan_delta=scan_delta,
                   darn_grid_stuff=darn_grid_stuff,
                   med_filter=med_filter,
                   model=model,
                   cache_geometry=cache_geometry)


def _freeze_model_arrays(model, min_bytes=1 << 20):
    """
    Mark the Emodel's large matrices read-only, so the workers share rather than copy them.
    The inversion only rebinds these attributes, so nothing should be writing through them.
    If something does, it raises instead of quietly diverging in one worker.
    Set FBI_FREEZE_MODEL=0 to skip.
    :param model: lompe Emodel
    :param min_bytes: smallest array worth freezing
    """

    if os.environ.get('FBI_FREEZE_MODEL') == '0':
        return

    for value in vars(model).values():
        if isinstance(value, np.ndarray) and value.nbytes >= min_bytes:
            try:
                value.flags.writeable = False
            except ValueError:
                # A view whose base is already read-only
                pass


def _reset_caches():
    """
    Clear the process-wide state set up by process().
    None of it is keyed on the grid or the apex epoch, so it can't be reused between calls.
    """

    global _worker_apex
    _worker_apex = None
    _gate_position_tables.clear()
    _shared.clear()
    readwrite.reset_geometry_cache()


def _lompe_one_scan(index):
    """
    Create the lompe fit for a single scan. This is what each worker runs.
    Everything but the index comes from _shared, which the workers inherit.
    :param index: int - position of this scan in range_times
    :return: lompe_extract() output, or None if no fit was made
    """

    scan_time = _shared['range_times'][index]
    all_data = _shared['all_data_iterable'][index]
    model = _shared['model']

    # Get the data in a format that Lompe likes
    sd_data, rids = prepare_lompe_inputs(_worker_apex, all_data, scan_time,
                                         _shared['scan_delta'], _shared['med_filter'])

    if sd_data is None:
        return None

    # The model is reused between scans, so drop the previous scan's data first
    model.clear_model()

    # Run lompe
    try:
        scan_lompe = run_lompe_model(sd_data, model)
    except IndexError:
        return None

    if scan_lompe is None:
        return None

    return lompe_extract(scan_lompe, _worker_apex, scan_time, _shared['darn_grid_stuff'],
                         rids, use_cache=_shared['cache_geometry'])


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
    # posterior=False skips Cmpost and Rmatrix, which cost ~4N^3 (N = grid_E.size) and are
    # never read by lompe_extract - only the model vector is used.
    # Needs the lompe fork at github.com/billetd/lompe. Stock lompe passes posterior
    # through to scipy.linalg.lstsq, which raises TypeError on every scan.
    try:
        model.run_inversion(l1=10, l2=0.1, lapack_driver='gelsy', posterior=False)
    except TypeError as err:
        raise TypeError(
            'run_inversion() rejected its arguments. FBI needs the lompe fork at '
            'github.com/billetd/lompe, which takes posterior= and provides '
            'prepare_biggrid().'
        ) from err
    return model