"""
This module contains code for creating lompe fits between a given timerange, using the already read in
SuperDARN data
"""
import apexpy
import lompe
import gc
import os
import time
import numpy as np
import FBI.grid as grid
import FBI.inversion as inversion
import FBI.los as los
import FBI.readwrite as readwrite
from FBI.fitacf import get_scan_times_widebeam, get_scan_times_old, record_time
from FBI.parallel import resolve_cores, forked_pool, bounded_imap, report_progress
from FBI.readwrite import lompe_extract, FBIWriter

# Everything the workers read. Filled in before forking, so the workers inherit it and a
# task is just a scan index.
_shared = {}


def process(all_data, timerange, lompe_dir, cores=1, med_filter=True, scandelta_override=None, range_times=None):
    """
    :param all_data: list[list[dict]] - Records of each radar, read in with fitacf.read_fitacfs()
    :param timerange: list[datetime] - Start and end times
    :param lompe_dir: str - Directory to save FBI output file
    :param cores: int - Number of worker processes. None uses every CPU available. Choose 1
                  to fit in this process with no pool, for profiling or debugging.
    :param med_filter: True or False - Median filter the data before putting into Lompe
    :param scandelta_override: int - Time in seconds to gather data around scans
    :param range_times: list[datetime] - Custom "scan" intervals
    """

    cores = resolve_cores(cores)

    # Clear anything left by a previous call, e.g. from extras.process_date()
    _shared.clear()

    # If no "range_times" is given, work it out based on input data
    # May produce silly scans if there is a mix of normal scan and other modes. Be cautious and only use
    # if you know what data is going in.
    if not range_times:
        # Get scan times within timerange, based on whichever radar started earlier
        # Run the old way if using "scanning data", or the new way if using widebeam data
        if all_data[0][1]['scan'] == 0:
            _, range_times, scan_delta = get_scan_times_old(all_data, timerange)
        else:
            range_times, scan_delta = get_scan_times_widebeam(all_data, timerange)

    # Override scan_delta here to integreate more data per scan
    if scandelta_override is not None:
        scan_delta = scandelta_override

    n_total = len(range_times)
    cores = min(cores, max(1, n_total))  # No point in more workers than scans

    # Initialise an apexpy object, for magnetic transforms
    apex = apexpy.Apex(range_times[0], refh=300)

    # Get the Superdarn grid
    # This is used for plotting purposes later, e.g. shading darker vectors where there is data
    darn_grid_stuff = grid.sdarn_grid(apex)

    # Cut down lompe model on a grid encompassing the SuperDARN Canada PolarDARNs, with the
    # regularisation used for every fit
    model = inversion.Model(grid.lompe_grid_canada(apex), l1=10, l2=0.1, ew_regularization_limit=(50, 75),
                            threads=cores)
    del apex  # No longer needed

    # Build everything the workers need before forking, so they share it rather than
    # each being sent a copy
    _prime_shared_state(all_data, range_times, scan_delta, darn_grid_stuff, med_filter, model)

    del darn_grid_stuff, model
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
        # Drop the model and data before the caller moves on to the next chunk
        _shared.clear()
        gc.collect()


def _prime_shared_state(all_data, range_times, scan_delta, darn_grid_stuff, med_filter, model):
    """
    Set up the state the workers read, before the pool is forked
    :param all_data: list[list[dict]] - Records of each radar
    :param range_times: list[datetime] - scan times
    :param scan_delta: int - seconds of data to gather around each scan
    :param darn_grid_stuff: dict from FBI.grid.sdarn_grid()
    :param med_filter: True or False
    :param model: FBI.inversion.Model
    """

    # Must come after the model is built. The model makes its own Apex at epoch 2015, and
    # apexpy holds the epoch in Fortran state shared by every Apex in the process, so
    # building ours last puts the epoch back to the one the data wants.
    apex = apexpy.Apex(range_times[0], refh=300)

    # Everything saved for a scan that doesn't depend on the data
    output_geometry = readwrite.output_geometry(model, apex, darn_grid_stuff)

    # Every gate any scan will use
    times = [[record_time(record) for record in radar] for radar in all_data]
    gates = los.Gates()
    for scan_time in range_times:
        for table in los.scan_tables(los.scan_windows(all_data, times, scan_time, scan_delta),
                                     scan_time, scan_delta):
            gates.table(*table)

    # Everything about those gates the fits need, so a scan only has to look it up
    gate_geometry = _gate_geometry(gates, model, apex)

    _freeze_arrays([*vars(model).values(), *gate_geometry.values(), *output_geometry.values(),
                    *output_geometry['v_matrix_model'], *output_geometry['v_matrix_darngrid']])

    _shared.update(all_data=all_data,
                   times=times,
                   range_times=range_times,
                   scan_delta=scan_delta,
                   med_filter=med_filter,
                   model=model,
                   gates=gates,
                   gate_geometry=gate_geometry,
                   output_geometry=output_geometry)


def _gate_geometry(gates, model, apex):
    """
    The parts of a fit that only depend on which gate a point is from
    :param gates: FBI.los.Gates
    :param model: FBI.inversion.Model
    :param apex: apexpy.Apex object
    :return: dict. 'row' is each gate's row in the others, or -1 if it is outside the area data
             are fitted in.
    """

    in_fit = model.biggrid.ingrid(gates.lon, gates.lat)
    row = np.full(gates.size, -1)
    row[in_fit] = np.arange(np.count_nonzero(in_fit))
    lon, lat = gates.lon[in_fit], gates.lat[in_fit]
    f1, f2 = apex.basevectors_qd(lat, lon, 300, coords='geo')
    mlat, mlon = apex.geo2apex(lat, lon, 300)

    return {'row': row,
            # Velocity along the line of sight from the model vector
            'G': model.los_matrix(lon, lat, gates.le[in_fit], gates.ln[in_fit]),
            'mlat': mlat, 'mlon': mlon, 'f1': f1, 'f2': f2}


def _freeze_arrays(values, min_bytes=1 << 20):
    """
    Mark large matrices read-only, so the workers share rather than copy them.
    Nothing should be writing through them. If something does, it raises instead of
    quietly diverging in one worker.
    Set FBI_FREEZE_MODEL=0 to skip.
    :param values: objects to freeze, anything that isn't an array is skipped
    :param min_bytes: smallest array worth freezing
    """

    if os.environ.get('FBI_FREEZE_MODEL') == '0':
        return

    for value in values:
        if isinstance(value, np.ndarray) and value.nbytes >= min_bytes:
            try:
                value.flags.writeable = False
            except ValueError:
                # A view whose base is already read-only
                pass


def _lompe_one_scan(index):
    """
    Create the lompe fit for a single scan. This is what each worker runs.
    Everything but the index comes from _shared, which the workers inherit.
    :param index: int - position of this scan in range_times
    :return: lompe_extract() output, or None if no fit was made
    """

    scan_time = _shared['range_times'][index]
    scan_delta = _shared['scan_delta']
    model, gates, geometry = _shared['model'], _shared['gates'], _shared['gate_geometry']

    # Line-of-sight data
    windows = los.scan_windows(_shared['all_data'], _shared['times'], scan_time, scan_delta)
    gate, vlos, vlos_err, rid = los.scan_los(gates, windows, scan_time, scan_delta, med_filter=_shared['med_filter'])

    # Same data selection as lompe: drop NaNs, then anything outside biggrid
    row = geometry['row'][gate]
    used = np.isfinite(vlos) & (row >= 0)
    if used.sum() <= 1:
        return None
    gate, row, vlos, vlos_err, rid = gate[used], row[used], vlos[used], vlos_err[used], rid[used]
    lon, lat = gates.lon[gate], gates.lat[gate]

    # Run lompe
    m = model.solve(geometry['G'], row, lon, lat, vlos, vlos_err)

    los_data = {'v_e_geo': vlos * gates.le[gate], 'v_n_geo': vlos * gates.ln[gate],
                'mlats': geometry['mlat'][row], 'mlons': geometry['mlon'][row],
                'f': (geometry['f1'][:, row], geometry['f2'][:, row]), 'rids': rid}

    return lompe_extract(_shared['output_geometry'], m, los_data, scan_time)


def prepare_lompe_inputs(all_data, scan_time, scan_delta, med_filter):
    """
    Make a lompe Data object of the SuperDARN line-of-sight velocities of a scan
    :param all_data: list[list[dict]] - Records of each radar, from fitacf.read_fitacfs()
    :param scan_time: datetime - Time of the scan
    :param scan_delta: float - Time in seconds to gather data around the scan
    :param med_filter: True or False - Median filter the data
    :return: (lompe.Data, station id of each point), or (None, None) if there is no data
    """

    gates = los.Gates()
    gate, vlos, vlos_err, rid = los.scan_los(gates, los.radar_windows(all_data), scan_time, scan_delta,
                                             med_filter=med_filter)

    # lompe wants the speed, with the direction in the LOS unit vector
    sign = np.sign(vlos)
    coords = np.vstack((gates.lon[gate], gates.lat[gate]))
    los_vectors = np.vstack((sign * gates.le[gate], sign * gates.ln[gate]))

    # Make the Lompe data object
    try:
        sd_data = lompe.Data(np.abs(vlos), coordinates=coords, LOS=los_vectors,
                             datatype='convection', error=vlos_err, iweight=1.0)
    except AttributeError:
        print('No data in this scan for some reason. Skipping...')
        return None, None

    return sd_data, rid
