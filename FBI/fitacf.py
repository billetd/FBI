"""
This module contains code for reading in, and handling, SuperDARN fitacf files
"""

import pydarn
import gc
import datetime as dt
import numpy as np
import math
from FBI.parallel import resolve_cores, forked_pool
from FBI.utils import find_indexes_within_time_range


def sdarnreadmulti(fitacf_file, start=None, end=None):
    """
    Read one fitacf file, keeping only the records and keys lompe needs
    :param fitacf_file: str, path to the file
    :param start: Start time to store, default the start of the file
    :param end: End time to store, default the end of the file
    :return: list[dict], one dictionary per record
    """

    print('Reading: ' + fitacf_file)

    # sdarnread = pydarn.SuperDARNRead(fitacf_file)
    # fitacf_data = sdarnread.read_fitacf()
    fitacf_data, _ = pydarn.read_fitacf(fitacf_file)

    # Keep only keys which are required for lompe
    keys_to_keep = ['time.yr', 'time.mo', 'time.dy', 'time.hr', 'time.mt', 'time.sc',
                    'time.us', 'scan', 'bmnum', 'stid', 'slist', 'gflg', 'rsep', 'frang',
                    'v', 'v_e', 'nrang']
    new_fitacf_data = []
    for d in fitacf_data:
        new_dict = {}
        for key in keys_to_keep:
            new_dict[key] = d.get(key)  # This will return None if key is missing
        new_fitacf_data.append(new_dict)
    del fitacf_data

    # Get all times in the file
    rid_record_times = [dt.datetime(new_fitacf_data[x]['time.yr'], new_fitacf_data[x]['time.mo'],
                                    new_fitacf_data[x]['time.dy'], new_fitacf_data[x]['time.hr'],
                                    new_fitacf_data[x]['time.mt'], new_fitacf_data[x]['time.sc'],
                                    new_fitacf_data[x]['time.us'])
                        for x in range(0, len(new_fitacf_data))]

    if start is None:
        start = rid_record_times[0]
    if end is None:
        end = rid_record_times[-1]

    filtered_indices = [i for i, time in enumerate(rid_record_times) if start <= time <= end]
    new_new_fitacf_data = [new_fitacf_data[index] for index in filtered_indices]

    gc.collect()
    return new_new_fitacf_data


def read_fitacfs(fitacf_files, cores=None, start=None, end=None):
    """
    Reads fitacf files into one big list [number of files] of list [number of records]
    of dictionaries [fitacf variables]
    :param fitacf_files: list[str], list of fitacf files to read. Must be a list/array, even if just one file
    :param cores: int, number of worker processes. None uses every CPU available. 1 reads
                  in this process, with no pool.
    :param start: Start time to store, default the start of the file
    :param end: End time to store, default the end of the file
    :return:
    """

    cores = min(resolve_cores(cores), max(1, len(fitacf_files)))

    if cores == 1:
        all_data = [sdarnreadmulti(inp, start, end) for inp in fitacf_files]
    else:
        # One file per worker
        with forked_pool(cores) as pool:
            all_data = pool.starmap(sdarnreadmulti,
                                    [(inp, start, end) for inp in fitacf_files])

    all_data = [x for x in all_data if x and x is not None]  # Gets rid of empty lists and None's
    return all_data


def all_data_make_iterable(all_data, range_times, scan_delta):
    """
    With the data read in from read_fitacfs(), remove the data that doesn't fall within our timerange,
    and turn it into an iterable
    :param all_data:
    :param range_times:
    :param scan_delta:
    :return:
    """

    all_data_iterable = []

    # Get the times of all records for each radar
    record_times = []
    for file_index in range(len(all_data)):
        rid_record_times = [dt.datetime(all_data[file_index][x]['time.yr'], all_data[file_index][x]['time.mo'],
                                        all_data[file_index][x]['time.dy'], all_data[file_index][x]['time.hr'],
                                        all_data[file_index][x]['time.mt'], all_data[file_index][x]['time.sc'],
                                        all_data[file_index][x]['time.us'])
                            for x in range(0, len(all_data[file_index]))]
        record_times.append(rid_record_times)

    for counter, scan_time in enumerate(range_times):

        all_radars_this_scan = []
        for file_index in range(len(all_data)):

            # Get the indexes for the records which are within scan_delta of scan_time
            # scan_delta multiplied by three to account for median filtering later (needs scan before and after)
            records_in_scan = find_indexes_within_time_range(record_times[file_index], scan_time,
                                                             catchtime=((scan_delta * 3) / 2))  # protect median filter
            these_records = [all_data[file_index][record] for record in records_in_scan]
            if these_records:  # prevent occurences where no records are found creeping in
                all_radars_this_scan.append(these_records)

        all_data_iterable.append(all_radars_this_scan)

    return all_data_iterable


def get_scan_times_widebeam(all_data, timerange):

    all_times = [
        [dt.datetime(record["time.yr"],
                     record["time.mo"],
                     record["time.dy"],
                     record["time.hr"],
                     record["time.mt"],
                     record["time.sc"],
                     record["time.us"])
         for record in all_data[radar]
         ]
        for radar in range(len(all_data))
    ]

    unique_times = [set(radar) for radar in all_times]
    sorted_times = [sorted(radar) for radar in unique_times]

    # Scan lengths
    diffs = [[(t2 - t1) for t1, t2 in zip(radar, radar[1:])] for radar in sorted_times]

    # Remove empty lists (files with incomplete records, or single scanes
    diffs = [sublist for sublist in diffs if sublist]

    # Check for data gaps, which would manifest as a diff much bigger than all the others
    # The mimimum diff should be close-ish to the average. Use that as a first "best guess"
    # Any diff that is off from 50% of the min is considered a delayed scan, and is removed when considering the avg
    # Need to subtract min from diffs and check to see which are over threshold
    min_diffs = [min(radar) for radar in diffs]
    diffs_threshold = [0.5*d for d in min_diffs]
    filtered_diffs = [[d for d in diffs[i] if d - min_diffs[i] < diffs_threshold[i]]
                      for i in range(len(diffs_threshold))]

    # Get average scan lengths
    average_diffs = [sum(radar, radar[0]-radar[0]) / len(radar) for radar in filtered_diffs]
    scan_delta = sum(average_diffs, average_diffs[0]-average_diffs[0]) / len(average_diffs)

    # Should start with the first record after the start of the timerange given
    # This avoids empty records
    first_time = min(
        (dt for sublist in unique_times for dt in sublist if dt > timerange[0]),
        default=None
    )
    range_times = [first_time + i * scan_delta
                   for i in range(int((timerange[1] - first_time) / scan_delta) + 1)]

    return range_times, scan_delta


def get_scan_times_old(all_data, timerange):
    """

    :param all_data: list[list[dict]], list of list of dictionaries containing all the fitacf data read in from
    read_fitacfs()
    :param timerange: list[datetime], list of two datetime objects denoting the start and end times
    :return: scan_times, list[datetime], times where the scans start, based on whichever radar starts first
    :return: range_times, list[datetime], times of scans that fall within timerange
    :return: range_times, list[datetime], times of scans that fall within timerange
    :return: scan_delta, float, time difference between successive scans
    """

    # Determine which radar starts the soonest
    # Will use that radars times as a baseline for all our scans
    start_times = [
        dt.datetime(data[0]['time.yr'], data[0]['time.mo'], data[0]['time.dy'],
                    data[0]['time.hr'], data[0]['time.mt'], data[0]['time.sc'],
                    data[0]['time.us'])
        for data in all_data
    ]

    # Find the index of the earliest datetime in the list
    earliest_radar = start_times.index(min(start_times))

    # Indexes where scan flag is 1
    scan_indexes = [index for index, data_dict in enumerate(all_data[earliest_radar]) if
                    data_dict.get("scan") == 1]

    # If no scans found, then use beam number
    # This might break weird scanning modes
    if not scan_indexes:
        scan_indexes = [index for index, data_dict in enumerate(all_data[earliest_radar]) if
                        data_dict.get("bmnum") == 0]

    # Special case for the imaging mode data
    # Only get the unique times if there as many records as scan,
    # Else only get the times where the scan flag is 1
    if len(all_data[earliest_radar]) == len(scan_indexes):
        scan_times = sorted(set(
            dt.datetime(entry["time.yr"], entry["time.mo"], entry["time.dy"], entry["time.hr"], entry["time.mt"],
                        entry["time.sc"], entry["time.us"]
                        )
            for entry in all_data[earliest_radar]
        ))
    else:
        scan_times = [
            dt.datetime(all_data[earliest_radar][index]["time.yr"], all_data[earliest_radar][index]["time.mo"],
                        all_data[earliest_radar][index]["time.dy"], all_data[earliest_radar][index]["time.hr"],
                        all_data[earliest_radar][index]["time.mt"], all_data[earliest_radar][index]["time.sc"],
                        all_data[earliest_radar][index]["time.us"]
                        )
            for index in scan_indexes]

    # Restrict to times that fall within timerange
    range_times = [time for time in scan_times if timerange[0] <= time <= timerange[1]]

    # Time difference between scans
    scan_delta = (scan_times[1]-scan_times[0]).seconds

    return scan_times, range_times, scan_delta


def median_filter_record(weighting_array, fitacf_data, record, max_beams):
    """
    Vectorised equivalent of median_filter(), evaluated for every range gate of a record
    in one pass.

    For a fixed record the 3x3 neighbourhood of (scan, beam) records is the same for every
    gate - only the 3-gate window moves. So the neighbourhood is unpacked once into dense
    arrays indexed by gate, and the filter becomes a 3-tap window slid along the gate axis.

    :param weighting_array: (3, 3, 3) array, indexed [scan_counter, gate_offset, beam_counter]
    :param fitacf_data: list[dict] of records for one radar
    :param record: index of the record to filter
    :param max_beams: highest bmnum in the file, plus one
    :return: (medians, passed), both of length nrang and indexed by range gate.
             medians is the filtered velocity (NaN where no scatter was found), passed is
             True where the weighting score was beaten, i.e. where median_filter() would
             have returned a value rather than [].
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

            beam_diff = beam_counter - 1  # This allows us to get the correct records
            index = scan + beam_diff

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


def median_filter(weighting_array, fitacf_data, record, max_beams, gate):
    """

    :param weighting_array:
    :param fitacf_data:
    :param record:
    :param max_beams:
    :param gate:
    :return:
    """

    # Score to beat when summing scatter in range gates. Will be halved if on a beam/range edge.
    weight_score = 24
    # weight_score = 6

    # Total number of records in this file
    n_recs = len(fitacf_data)

    # Total number of range gates. Assumes constant in file (might break?)
    max_range = fitacf_data[record]['nrang']

    # Previous and next scan
    scans = [record - max_beams, record, record + max_beams]

    # Current, left, and right beams
    beams = np.array([fitacf_data[record]['bmnum'] - 1,
                      fitacf_data[record]['bmnum'],
                      fitacf_data[record]['bmnum'] + 1])
    if np.logical_or(beams[0] < 0, beams[2] > max_beams):
        weight_score -= 3  # Reduce weight score by number of lost cells

    # Current, up, and down ranges
    gates = np.array([gate - 1, gate, gate + 1])
    if np.logical_or(gates[0] < 0, gates[2] > max_range):
        weight_score -= 3  # Reduce weight score by number of lost cells

    cumulative_weight = 0
    vels = []

    # Iterate over previous, current, and next scans
    for scan_counter, scan in enumerate(scans):
        if np.logical_and(scan >= 0, scan < n_recs):  # Check not before or after start/end of records
            scatter = np.zeros([3, 3])

            # Iterating over left, current, rigt beams
            for beam_counter, beam in enumerate(beams):
                if np.logical_and(beam >= 0, beam < max_beams):

                    beam_diff = beam_counter - 1  # This allows us to get the correct records

                    # Have to check again if there is actually any data
                    try:
                        current_beam_slist = fitacf_data[scan + beam_diff]['slist']
                    except (KeyError, IndexError) as err:
                        continue
                    current_beam_gscat = fitacf_data[scan + beam_diff]['gflg']

                    # print(scan_counter, beam_counter)
                    # print(fitacf_data[scan + beam_diff]['bmnum'])

                    # Check the gates are in slist
                    isin = np.isin([gates[0], gates[1], gates[2]], current_beam_slist)
                    isin_indexes = np.where(isin)[0]

                    # Check the positions are not ground scatter
                    if isin_indexes.size != 0:
                        slist_indexes = np.where(np.isin(current_beam_slist, gates))[0]
                        gscat = np.where(current_beam_gscat[slist_indexes] == 0)
                        scatter[isin_indexes[gscat], beam_counter] = 1  # Indexing the current beam
                        vels.extend(fitacf_data[scan + beam_diff]['v'][slist_indexes[gscat]])

            # Sum the weights and add it to the counter
            cumulative_weight += np.sum(scatter * weighting_array[scan_counter])

    # Finally, median filter if we beat the weighting score, otherwise return []
    if cumulative_weight > weight_score:
        return np.median(vels)
    else:
        return []


def fitacf_get_k_vector_circle(radlat, radlon, radmlat, radmlon,
                              lats, lons, mlats, mlons, vlos):
    """

    :param radlat:  float  - radar geographic latitude
    :param radlon:  float  - radar geographic longitude
    :param radmlat: float  - radar magnetic latitude
    :param radmlon: float  - radar magnetic longitude
    :param lats:    ndarray - gate geographic latitudes
    :param lons:    ndarray - gate geographic longitudes
    :param mlats:   ndarray - gate magnetic latitudes
    :param mlons:   ndarray - gate magnetic longitudes
    :param vlos:    ndarray - signed line-of-sight velocities (m/s)
    :return: le, ln, le_mag, ln_mag, ve_geo, vn_geo, ve_mag, vn_mag — all ndarray
    """

    # Graciously adapted from Evan's code (invmag.pro), adding vectorisation

    # Geographic
    aside = 90.0 - radlat
    cos_a = math.cos(math.radians(aside))
    sin_a = math.sin(math.radians(aside))
    cside = 90.0 - lats
    Bangle = radlon - lons

    # Haversine formula
    arg = (cos_a * np.cos(np.radians(cside))
           + sin_a * np.sin(np.radians(cside)) * np.cos(np.radians(Bangle)))
    arg = np.clip(arg, -1.0, 1.0)  # guard against float rounding outside [-1, 1]
    bside = np.degrees(np.arccos(arg))

    numer = cos_a - np.cos(np.radians(bside)) * np.cos(np.radians(cside))
    denom = np.sin(np.radians(bside)) * np.sin(np.radians(cside))
    denom = np.where(np.abs(denom) < 1e-10, 1e-10, denom)  # guard against gate-at-radar
    arg2 = np.clip(numer / denom, -1.0, 1.0)
    Aangle_geo = np.degrees(np.arccos(arg2))
    az_geo = np.where(Bangle < 0, -Aangle_geo, Aangle_geo)
    az_geo = np.where(np.isnan(az_geo), 0.0, az_geo)

    le = np.sign(vlos) * np.sin(np.radians(az_geo))
    ln = np.sign(vlos) * np.cos(np.radians(az_geo))
    ve_geo = vlos * np.sin(np.radians(az_geo))
    vn_geo = vlos * np.cos(np.radians(az_geo))

    # Same but in magnetic
    aside_m = 90.0 - radmlat
    cos_am = math.cos(math.radians(aside_m))
    sin_am = math.sin(math.radians(aside_m))
    cside_m = 90.0 - mlats
    Bangle_m = radmlon - mlons

    # Haversine
    arg_m = (cos_am * np.cos(np.radians(cside_m))
             + sin_am * np.sin(np.radians(cside_m)) * np.cos(np.radians(Bangle_m)))
    arg_m = np.clip(arg_m, -1.0, 1.0)
    bside_m = np.degrees(np.arccos(arg_m))

    numer_m = cos_am - np.cos(np.radians(bside_m)) * np.cos(np.radians(cside_m))
    denom_m = np.sin(np.radians(bside_m)) * np.sin(np.radians(cside_m))
    denom_m = np.where(np.abs(denom_m) < 1e-10, 1e-10, denom_m)
    arg2_m = np.clip(numer_m / denom_m, -1.0, 1.0)
    Aangle_mag = np.degrees(np.arccos(arg2_m))
    az_mag = np.where((Bangle_m < 0) & (np.abs(Bangle_m) < 180), -Aangle_mag, Aangle_mag)
    az_mag = np.where(np.isnan(az_mag), 0.0, az_mag) # Check for cases when az = NAN rather than zero

    le_mag = np.sign(vlos) * np.sin(np.radians(az_mag))
    ln_mag = np.sign(vlos) * np.cos(np.radians(az_mag))
    ve_mag = vlos * np.sin(np.radians(az_mag))
    vn_mag = vlos * np.cos(np.radians(az_mag))

    return le, ln, le_mag, ln_mag, ve_geo, vn_geo, ve_mag, vn_mag


