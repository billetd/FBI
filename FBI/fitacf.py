"""
This module contains code for reading in, and handling, SuperDARN fitacf files
"""

import pydarn
import gc
import os
import datetime as dt
import numpy as np
import math
from FBI.utils import find_indexes_within_time_range


def read_fitacfs(fitacf_files, cores=1):
    """
    Reads fitacf files into one big list [number of files] of list [number of records]
    of dictionaries [fitacf variables]
    :param fitacf_files: list[str], list of fitacf files to read
    :param cores: int, number of cores to use when parallel reading files. Defaults to zero.
    :return:
    """

    os.environ['RAY_DEDUP_LOGS'] = '0'
    import ray

    ray.init(num_cpus=cores)

    @ray.remote
    def sdarnreadmulti(fitacf_file):
        print('Reading: ' + fitacf_file)
        sdarnread = pydarn.SuperDARNRead(fitacf_file)
        fitacf_data = sdarnread.read_fitacf()

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
        gc.collect()
        return new_fitacf_data

    # Read in all the files
    all_data = ray.get([sdarnreadmulti.remote(inp) for inp in fitacf_files])
    ray.shutdown()
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
    for counter, scan_time in enumerate(range_times):
        all_radars_this_scan = []
        for file_index in range(len(all_data)):

            # Get the indexes for the records which are within scan_delta of scan_time
            # scan_delta multiplied by three to account for median filtering later (needs scan before and after)
            record_times = [dt.datetime(all_data[file_index][x]['time.yr'], all_data[file_index][x]['time.mo'],
                                        all_data[file_index][x]['time.dy'], all_data[file_index][x]['time.hr'],
                                        all_data[file_index][x]['time.mt'], all_data[file_index][x]['time.sc'],
                                        all_data[file_index][x]['time.us'])
                            for x in range(0, len(all_data[file_index]))]
            records_in_scan = find_indexes_within_time_range(record_times, scan_time,
                                                             catchtime=((scan_delta * 3) / 2))  # protect median filter
            these_records = [all_data[file_index][record] for record in records_in_scan]
            if these_records:  # prevent occurences where no records are found creeping in
                all_radars_this_scan.append(these_records)

            if file_index == 3:
                pass

        all_data_iterable.append(all_radars_this_scan)

    return all_data_iterable


def get_scan_times(all_data, timerange):
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

    # Get scan flags
    scan_flags = [entry["scan"] for entry in all_data[earliest_radar]]
    beams = [entry["bmnum"] for entry in all_data[earliest_radar]]

    # Indexes where scan flag is 1
    scan_indexes = [index for index, data_dict in enumerate(all_data[earliest_radar]) if
                    data_dict.get("scan") == 1]

    # If no scans found, then use beam number
    # This might break weird scanning modes
    if not scan_indexes:
        scan_indexes = [index for index, data_dict in enumerate(all_data[earliest_radar]) if
                        data_dict.get("bmnum") == 0]

    # All the times from all records
    all_times = sorted(
        dt.datetime(entry["time.yr"], entry["time.mo"], entry["time.dy"], entry["time.hr"], entry["time.mt"],
                    entry["time.sc"], entry["time.us"]
                    )
        for entry in all_data[earliest_radar]
    )

    # Special case for the imaging mode data
    # Only get the unique times if there as many records as scan,
    # Else only get the times where the scan flag is 1
    if len(all_times) == len(scan_indexes):
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


def median_filter(fitacf_data, record, max_beams, gate):
    """

    :param fitacf_data:
    :param record:
    :param max_beams:
    :param gate:
    :return:
    """

    # Score to beat when summing scatter in range gates. Will be halved if on a beam/range edge.
    weight_score = 24

    weighting_array = np.array([[[1, 1, 1], [1, 2, 1], [1, 1, 1]],
                                [[2, 2, 2], [2, 4, 2], [2, 2, 2]],
                                [[1, 1, 1], [1, 2, 1], [1, 1, 1]]
                                ])

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
                    except KeyError:
                        continue
                    current_beam_gscat = fitacf_data[scan + beam_diff]['gflg']

                    # Check the gates are in slist
                    isin = np.isin([gates[0], gates[1], gates[2]], current_beam_slist)
                    isin_indexes = np.where(isin)[0]

                    # Check the positions are not ground scatter
                    if isin_indexes.size != 0:
                        slist_indexes = np.where(np.isin(current_beam_slist, gates))[0]
                        gscat = np.where(current_beam_gscat[slist_indexes] == 0)
                        scatter[isin_indexes[gscat], beam_counter] = 1  # Indexing the current beam
                        vels = np.append(vels, fitacf_data[scan + beam_diff]['v'][slist_indexes[gscat]])

            # Sum the weights and add it to the counter
            cumulative_weight += np.sum(scatter * weighting_array[scan_counter])

    # Finally, median filter if we beat the weighting score, otherwise return []
    if cumulative_weight > weight_score:
        return np.median(vels)
    else:
        return []


def fitacf_get_k_vector_circle(apex, stid, lat, lon, mlat, mlon, v_los):
    """

    :param apex:
    :param stid:
    :param lat:
    :param lon:
    :param mlat:
    :param mlon:
    :param v_los:
    :return:
    """

    # Get position of radar in geographic from hdw files in pyDARN
    radlat = pydarn.SuperDARNRadars.radars[stid].hardware_info.geographic.lat
    radlon = pydarn.SuperDARNRadars.radars[stid].hardware_info.geographic.lon
    radmlat, radmlon = apex.geo2apex(radlat, radlon, 300)

    # Graciously copied from Evan's code (invmag.pro)
    # Geographic
    api = 4 * math.atan(1.0)
    aside = 90 - radlat
    cside = 90 - lat
    Bangle = radlon - lon

    # Haversine formula
    arg = (math.cos(aside * api / 180.0) * math.cos(cside * api / 180.0) +
           math.sin(aside * api / 180.0) * math.sin(cside * api / 180.0) * math.cos(Bangle * api / 180.0)
           )
    bside = math.acos(arg) * 180.0 / api

    arg2 = ((math.cos(aside * api / 180.0) - math.cos(bside * api / 180.0) * math.cos(cside * api / 180.0)) /
            (math.sin(bside * api / 180.0) * math.sin(cside * api / 180.0)))

    Aangle = math.acos(arg2) * 180.0 / api

    if Bangle < 0:
        Aangle = -Aangle

    az = Aangle

    # Check for cases when az = NAN rather than zero
    if math.isnan(az) is True:
        az = 0

    # v_e = v_los * np.sin(np.radians(az))  # Negative because positive towards radar
    # v_n = v_los * np.cos(np.radians(az))
    # kvect = np.degrees(math.atan2(v_e, v_n))

    ve_geo = v_los * np.sin(np.radians(az))
    vn_geo = v_los * np.cos(np.radians(az))
    le_current = np.sign(v_los) * np.sin(np.radians(az))
    ln_current = np.sign(v_los) * np.cos(np.radians(az))

    # Magnetic
    api = 4 * math.atan(1.0)
    aside = 90 - radmlat
    cside = 90 - mlat
    Bangle = radmlon - mlon

    # Haversine formula
    arg = (math.cos(aside * api / 180.0) * math.cos(cside * api / 180.0) +
           math.sin(aside * api / 180.0) * math.sin(cside * api / 180.0) * math.cos(Bangle * api / 180.0)
           )
    bside = math.acos(arg) * 180.0 / api

    arg2 = ((math.cos(aside * api / 180.0) - math.cos(bside * api / 180.0) * math.cos(cside * api / 180.0)) /
            (math.sin(bside * api / 180.0) * math.sin(cside * api / 180.0)))

    Aangle = math.acos(arg2) * 180.0 / api

    if Bangle < 0 and abs(Bangle) < 180:
        Aangle = -Aangle

    az = Aangle

    # Check for cases when az = NAN rather than zero
    if math.isnan(az) is True:
        az = 0

    ve_mag = v_los * np.sin(np.radians(az))
    vn_mag = v_los * np.cos(np.radians(az))
    le_mag_current = np.sign(v_los) * np.sin(np.radians(az))
    ln_mag_current = np.sign(v_los) * np.cos(np.radians(az))

    return le_current, ln_current, le_mag_current, ln_mag_current, ve_geo, vn_geo, ve_mag, vn_mag
