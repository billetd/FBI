"""
Line-of-sight velocities from fitacf records, and the geometry of the radar gates they come from
"""
import math
import numpy as np
import pydarn
from pydarn.utils.coordinates import gate2geographic_location
from FBI.fitacf import median_filter_record, record_time, time_window

# For median filtering
_WEIGHTING_ARRAY = np.array([[[1, 1, 1], [1, 2, 1], [1, 1, 1]],
                             [[2, 2, 2], [2, 4, 2], [2, 2, 2]],
                             [[1, 1, 1], [1, 2, 1], [1, 1, 1]]])


class Gates:
    """
    Position and line-of-sight direction of every radar gate looked up so far, numbered in the
    order they were added. Gates are added a table at a time, one table per
    (radar_id, max_beams, nrang, rsep, frang).
    """

    def __init__(self):
        self._tables = {}
        self.lat, self.lon = np.empty(0), np.empty(0)
        self.le, self.ln = np.empty(0), np.empty(0)

    @property
    def size(self):
        return self.lat.size

    def table(self, radar_id, max_beams, nrang, rsep, frang):
        """
        Numbers of the gates of a table, which is added if it isn't already
        :return: array of int, shaped (max_beams, nrang)
        """

        key = (radar_id, int(max_beams), int(nrang), int(rsep), int(frang))
        numbers = self._tables.get(key)
        if numbers is None:
            lat, lon, le, ln = _gate_table(*key)
            numbers = np.arange(self.size, self.size + lat.size).reshape(key[1], key[2])
            self._tables[key] = numbers
            self.lat, self.lon = np.concatenate((self.lat, lat)), np.concatenate((self.lon, lon))
            self.le, self.ln = np.concatenate((self.le, le)), np.concatenate((self.ln, ln))

        return numbers


def _gate_table(radar_id, max_beams, nrang, rsep, frang):
    """
    Geographic position and line-of-sight direction of every beam/gate of a radar. These have to
    be done in one call, as gate2geographic_location() iterates every gate until the slowest
    converges, so a gate's position depends on which others it is done with.
    :param radar_id: pydarn.RadarID
    :param max_beams: number of beams to tabulate
    :param nrang: number of range gates to tabulate
    :param rsep: range separation [km]
    :param frang: distance to the first range gate [km]
    :return: (lat, lon, le, ln) - Flattened over (max_beams, nrang). le and ln are the eastward and
             northward components of the unit vector from the gate towards the radar, the direction
             of a positive line-of-sight velocity.
    """

    beams, gates = np.meshgrid(np.arange(max_beams), np.arange(nrang), indexing='ij')
    lat, lon = gate2geographic_location(stid=radar_id, beam=beams.ravel(), range_gate=gates.ravel(),
                                        height=300, center=True, rsep=rsep, frang=frang)
    lat, lon = np.asarray(lat), np.asarray(lon)

    radar = pydarn.SuperDARNRadars.radars[radar_id].hardware_info.geographic
    azimuth = np.radians(los_azimuth(radar.lat, radar.lon, lat, lon))

    return lat, lon, np.sin(azimuth), np.cos(azimuth)


def los_azimuth(radlat, radlon, lats, lons):
    """
    Azimuth of the direction from each gate towards the radar
    Graciously adapted from Evan's code (invmag.pro), adding vectorisation
    :param radlat: float - Radar geographic latitude
    :param radlon: float - Radar geographic longitude
    :param lats: array - Gate geographic latitudes
    :param lons: array - Gate geographic longitudes
    :return: array - Azimuth [degrees east of north]
    """

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
    Aangle = np.degrees(np.arccos(arg2))
    azimuth = np.where(Bangle < 0, -Aangle, Aangle)

    return np.where(np.isnan(azimuth), 0.0, azimuth)


def scan_windows(all_data, times, scan_time, scan_delta):
    """
    Records of each radar around a scan, wide enough for the median filter to have the scans
    before and after
    :param all_data: list[list[dict]] - Records of each radar, from fitacf.read_fitacfs()
    :param times: list[list[datetime]] - Times of those records
    :param scan_time: datetime
    :param scan_delta: float - Seconds of data to gather around scans
    :return: list of (records, times), for radars with records around the scan
    """

    windows = []
    for records, radar_times in zip(all_data, times):
        window = time_window(radar_times, scan_time, scan_delta * 3 / 2)
        if window:
            windows.append((records[window.start:window.stop], radar_times[window.start:window.stop]))

    return windows


def radar_windows(all_data):
    """
    Treat all of each radar's records as one window
    :param all_data: list[list[dict]] - Records of each radar
    :return: list of (records, times)
    """

    return [(records, [record_time(record) for record in records]) for records in all_data]


def _scan_records(windows, scan_time, scan_delta):
    """
    The records of a scan that have scatter in them
    :param windows: list of (records, times) for each radar, from around the scan
    :param scan_time: datetime
    :param scan_delta: float - Seconds of data to gather around scans
    :return: generator of (records, index, max_beams, table) - The radar's records, the index of
             this record in them, and the Gates.table() arguments for its gates
    """

    for records, times in windows:
        radar_id = pydarn.RadarID(records[0]['stid'])
        max_beams = max([record['bmnum'] for record in records]) + 1

        for index in time_window(times, scan_time, scan_delta / 2):
            record = records[index]
            if record.get('slist') is None:
                continue
            table = (radar_id, max_beams, record['nrang'], record.get('rsep', 45), record.get('frang', 180))
            yield records, index, max_beams, table


def scan_tables(windows, scan_time, scan_delta):
    """
    The gate tables scan_los() needs for a scan, so they can be built in advance
    :return: set of Gates.table() arguments
    """

    return {table for _, _, _, table in _scan_records(windows, scan_time, scan_delta)}


def scan_los(gates, windows, scan_time, scan_delta, med_filter=False):
    """
    The line-of-sight velocities of one scan
    :param gates: Gates
    :param windows: list of (records, times) for each radar, from around the scan
    :param scan_time: datetime
    :param scan_delta: float - Seconds of data to gather around scans
    :param med_filter: True or False - Median filter the data
    :return: (gate, vlos, vlos_err, rid) - For each point, its gate number in gates, line-of-sight
             velocity (positive towards the radar), velocity error, and station id
    """

    gate, vlos, vlos_err, rid = [], [], [], []

    for records, index, max_beams, table in _scan_records(windows, scan_time, scan_delta):
        record = records[index]
        slist, v = record['slist'], record['v']
        numbers = gates.table(*table)[record['bmnum'], slist]

        # Only keep gates that are not ground scatter, have velocity below 2000m/s, and
        # are above range gate 10. This removes most erroneous data and near-range
        # (E-region) echos.
        keep = (record['gflg'] == 0) & (np.abs(v) <= 2000) & (slist > 10)

        # Median filtering. One pass covers every gate of the record.
        if med_filter:
            if not keep.any():
                continue
            medians, passed = median_filter_record(_WEIGHTING_ARRAY, records, index, max_beams)
            vel = medians[slist]
            # The median filter can fail if unreliable scatter is found
            keep &= passed[slist]
        else:
            vel = v

        # Velocities of exactly zero are dropped too
        keep &= vel != 0.0
        if not keep.any():
            continue

        gate.append(numbers[keep])
        vlos.append(vel[keep])
        vlos_err.append(np.abs(record['v_e'][keep]))
        rid.append(np.full(np.count_nonzero(keep), records[0]['stid']))

    if not gate:
        return np.empty(0, dtype=int), np.empty(0), np.empty(0), np.empty(0, dtype=int)

    return np.concatenate(gate), np.concatenate(vlos), np.concatenate(vlos_err), np.concatenate(rid)
