"""
This module contains code for reading in, and handling, SuperDARN fitacf files
"""

import bisect
import datetime as dt
import gc
import numpy as np
import pydarn
from FBI.parallel import resolve_cores, forked_pool

# Keys which are required for lompe
_KEYS = ['time.yr', 'time.mo', 'time.dy', 'time.hr', 'time.mt', 'time.sc', 'time.us',
         'scan', 'bmnum', 'stid', 'slist', 'gflg', 'rsep', 'frang', 'v', 'v_e', 'nrang']


def record_time(record):
    """
    :param record: dict - One fitacf record
    :return: datetime of the record
    """

    return dt.datetime(record['time.yr'], record['time.mo'], record['time.dy'], record['time.hr'],
                       record['time.mt'], record['time.sc'], record['time.us'])


def time_window(times, target_time, halfwidth):
    """
    Indexes of the times within halfwidth of target_time
    :param times: list[datetime] - Sorted times
    :param target_time: datetime
    :param halfwidth: float - [seconds]
    :return: range
    """

    return range(bisect.bisect_right(times, target_time - dt.timedelta(seconds=halfwidth)),
                 bisect.bisect_right(times, target_time + dt.timedelta(seconds=halfwidth)))


def sdarnreadmulti(fitacf_file, start=None, end=None):
    """
    Read one fitacf file, keeping only the records and keys lompe needs
    :param fitacf_file: str, path to the file
    :param start: Start time to store, default the start of the file
    :param end: End time to store, default the end of the file
    :return: list[dict], one dictionary per record
    """

    print('Reading: ' + fitacf_file)

    fitacf_data, _ = pydarn.read_fitacf(fitacf_file)
    records = [{key: record.get(key) for key in _KEYS} for record in fitacf_data]  # None if a key is missing
    del fitacf_data

    times = [record_time(record) for record in records]
    if start is None:
        start = times[0]
    if end is None:
        end = times[-1]
    records = [record for record, time in zip(records, times) if start <= time <= end]

    gc.collect()
    return records


def read_fitacfs(fitacf_files, cores=None, start=None, end=None):
    """
    Reads fitacf files into one big list [number of files] of list [number of records]
    of dictionaries [fitacf variables]
    :param fitacf_files: list[str], list of fitacf files to read. Must be a list/array, even if just one file
    :param cores: int, number of worker processes. None uses every CPU available. 1 reads
                  in this process, with no pool.
    :param start: Start time to store, default the start of the file
    :param end: End time to store, default the end of the file
    :return: list[list[dict]], without the files that had no records
    """

    cores = min(resolve_cores(cores), max(1, len(fitacf_files)))

    if cores == 1:
        all_data = [sdarnreadmulti(inp, start, end) for inp in fitacf_files]
    else:
        # One file per worker
        with forked_pool(cores) as pool:
            all_data = list(pool.map(sdarnreadmulti, fitacf_files, [start] * len(fitacf_files),
                                     [end] * len(fitacf_files)))

    return [x for x in all_data if x]


def get_scan_times_widebeam(all_data, timerange):
    """
    Scan times for widebeam data, which has no scan flag, spaced by the average time between integrations
    :param all_data: list[list[dict]] - From read_fitacfs()
    :param timerange: list[datetime] - Start and end times
    :return: range_times, list[datetime], times of scans that fall within timerange
    :return: scan_delta, timedelta, time between successive scans
    """

    unique_times = [set(record_time(record) for record in radar) for radar in all_data]
    sorted_times = [sorted(radar) for radar in unique_times]

    # Scan lengths
    diffs = [[(t2 - t1) for t1, t2 in zip(radar, radar[1:])] for radar in sorted_times]

    # Remove empty lists (files with incomplete records, or single scans)
    diffs = [sublist for sublist in diffs if sublist]

    # Check for data gaps, which would manifest as a diff much bigger than all the others
    # The mimimum diff should be close-ish to the average. Use that as a first "best guess"
    # Any diff that is off from 50% of the min is considered a delayed scan, and is removed when considering the avg
    min_diffs = [min(radar) for radar in diffs]
    diffs_threshold = [0.5 * d for d in min_diffs]
    filtered_diffs = [[d for d in diffs[i] if d - min_diffs[i] < diffs_threshold[i]]
                      for i in range(len(diffs_threshold))]

    # Get average scan lengths
    average_diffs = [sum(radar, radar[0] - radar[0]) / len(radar) for radar in filtered_diffs]
    scan_delta = sum(average_diffs, average_diffs[0] - average_diffs[0]) / len(average_diffs)

    # Should start with the first record after the start of the timerange given
    # This avoids empty records
    first_time = min((time for radar in unique_times for time in radar if time > timerange[0]), default=None)
    range_times = [first_time + i * scan_delta
                   for i in range(int((timerange[1] - first_time) / scan_delta) + 1)]

    return range_times, scan_delta


def get_scan_times_old(all_data, timerange):
    """
    Scan times for data with a scan flag, based on whichever radar starts first
    :param all_data: list[list[dict]] - From read_fitacfs()
    :param timerange: list[datetime] - Start and end times
    :return: scan_times, list[datetime], times where the scans start
    :return: range_times, list[datetime], times of scans that fall within timerange
    :return: scan_delta, int, seconds between successive scans
    """

    # Determine which radar starts the soonest
    # Will use that radars times as a baseline for all our scans
    start_times = [record_time(radar[0]) for radar in all_data]
    earliest = all_data[start_times.index(min(start_times))]

    # Indexes where scan flag is 1
    scan_indexes = [index for index, record in enumerate(earliest) if record.get("scan") == 1]

    # If no scans found, then use beam number
    # This might break weird scanning modes
    if not scan_indexes:
        scan_indexes = [index for index, record in enumerate(earliest) if record.get("bmnum") == 0]

    # Special case for the imaging mode data
    # Only get the unique times if there as many records as scan,
    # Else only get the times where the scan flag is 1
    if len(earliest) == len(scan_indexes):
        scan_times = sorted(set(record_time(record) for record in earliest))
    else:
        scan_times = [record_time(earliest[index]) for index in scan_indexes]

    # Restrict to times that fall within timerange
    range_times = [time for time in scan_times if timerange[0] <= time <= timerange[1]]

    # Time difference between scans
    scan_delta = (scan_times[1] - scan_times[0]).seconds

    return scan_times, range_times, scan_delta


def median_filter_record(weighting_array, fitacf_data, record, max_beams):
    """
    Median filter every range gate of a record, in one pass.

    For a fixed record the 3x3 neighbourhood of (scan, beam) records is the same for every
    gate - only the 3-gate window moves. So the neighbourhood is unpacked once into dense
    arrays indexed by gate, and the filter becomes a 3-tap window slid along the gate axis.

    :param weighting_array: (3, 3, 3) array, indexed [scan_counter, gate_offset, beam_counter]
    :param fitacf_data: list[dict] of records for one radar
    :param record: index of the record to filter
    :param max_beams: highest bmnum in fitacf_data, plus one
    :return: (medians, passed), both of length nrang and indexed by range gate.
             medians is the filtered velocity (NaN where no scatter was found), passed is
             True where the weighting score was beaten.
    """

    n_recs = len(fitacf_data)

    # Total number of range gates. Assumes constant in file (might break?)
    max_range = fitacf_data[record]['nrang']
    nrang = int(max_range)

    # Previous and next scan
    scans = [record - max_beams, record, record + max_beams]

    # Current, left, and right beams
    bmnum = fitacf_data[record]['bmnum']
    beams = [bmnum - 1, bmnum, bmnum + 1]

    # Score to beat when summing scatter in range gates. Reduced on beam/range edges.
    base_score = 24
    if beams[0] < 0 or beams[2] > max_beams:
        base_score -= 3  # Reduce weight score by number of lost cells

    # The range edge penalty depends on the gate, so it becomes an array
    gate_axis = np.arange(nrang)
    weight_score = np.full(nrang, base_score, dtype=np.int64)
    weight_score[(gate_axis - 1 < 0) | (gate_axis + 1 > max_range)] -= 3

    # Dense per-neighbour arrays over a gate axis padded by one on each side, so that
    # padded index p holds gate p - 1 and the gates -1 and nrang are addressable.
    npad = nrang + 2
    valid = np.zeros((3, 3, npad), dtype=bool)
    vel_pad = [[None] * 3 for _ in range(3)]
    vel_dtypes = []

    for scan_counter, scan in enumerate(scans):
        if not (0 <= scan < n_recs):  # Check not before or after start/end of records
            continue

        for beam_counter, beam in enumerate(beams):
            if not (0 <= beam < max_beams):
                continue

            index = scan + beam_counter - 1

            # Have to check again if there is actually any data
            try:
                nb_slist = fitacf_data[index]['slist']
            except (KeyError, IndexError):
                continue
            if nb_slist is None:
                continue

            nb_slist = np.asarray(nb_slist)
            if nb_slist.size == 0:
                continue

            nb_gflg = fitacf_data[index]['gflg']
            nb_v = np.asarray(fitacf_data[index]['v'])

            # Keep the non-ground-scatter gates that can fall inside a 3-gate window of
            # this record. Gates outside [-1, nrang] can never be reached.
            keep = (np.asarray(nb_gflg) == 0) & (nb_slist >= -1) & (nb_slist <= nrang)
            if not keep.any():
                continue

            positions = nb_slist[keep].astype(np.intp) + 1
            valid[scan_counter, beam_counter, positions] = True

            padded = np.zeros(npad, dtype=nb_v.dtype)
            padded[positions] = nb_v[keep]
            vel_pad[scan_counter][beam_counter] = padded
            vel_dtypes.append(nb_v.dtype)

    vel_dtype = np.result_type(*vel_dtypes) if vel_dtypes else np.dtype(np.float64)

    # Sum the weights of every occupied cell in the 3x3x3 neighbourhood of each gate, and
    # collect the velocities of those cells as columns for the median.
    cumulative_weight = np.zeros(nrang, dtype=np.float64)
    columns = []
    for scan_counter in range(3):
        weights = weighting_array[scan_counter]  # 3x3, [gate_offset, beam_counter]
        for beam_counter in range(3):
            padded = vel_pad[scan_counter][beam_counter]
            if padded is None:
                continue
            occupied = valid[scan_counter, beam_counter]
            for gate_offset in range(3):
                window = occupied[gate_offset:gate_offset + nrang]
                cumulative_weight += weights[gate_offset][beam_counter] * window
                columns.append((window, padded[gate_offset:gate_offset + nrang]))

    passed = cumulative_weight > weight_score

    # Per-gate median over however many cells were occupied. Filling the unused slots with
    # NaN lets a single sort handle all the different lengths at once, since NaNs sort last.
    medians = np.full(nrang, np.nan, dtype=vel_dtype)
    if columns:
        stack = np.full((nrang, len(columns)), np.nan, dtype=vel_dtype)
        counts = np.zeros(nrang, dtype=np.intp)
        for column, (window, values) in enumerate(columns):
            stack[window, column] = values[window]
            counts += window
        stack.sort(axis=1)
        rows = np.nonzero(counts > 0)[0]
        if rows.size:
            n = counts[rows]
            medians[rows] = 0.5 * (stack[rows, (n - 1) // 2] + stack[rows, n // 2])

    return medians, passed
