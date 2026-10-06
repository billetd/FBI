import numpy as np
import pytest
from FBI.fitacf import median_filter_record
from FBI.los import _WEIGHTING_ARRAY


def median_filter(weighting_array, fitacf_data, record, max_beams, gate):
    """The original one gate at a time median filter, which median_filter_record() vectorises"""

    weight_score = 24
    n_recs = len(fitacf_data)
    max_range = fitacf_data[record]['nrang']
    scans = [record - max_beams, record, record + max_beams]
    beams = np.array([fitacf_data[record]['bmnum'] - 1, fitacf_data[record]['bmnum'],
                      fitacf_data[record]['bmnum'] + 1])
    if np.logical_or(beams[0] < 0, beams[2] > max_beams):
        weight_score -= 3
    gates = np.array([gate - 1, gate, gate + 1])
    if np.logical_or(gates[0] < 0, gates[2] > max_range):
        weight_score -= 3

    cumulative_weight = 0
    vels = []
    for scan_counter, scan in enumerate(scans):
        if np.logical_and(scan >= 0, scan < n_recs):
            scatter = np.zeros([3, 3])
            for beam_counter, beam in enumerate(beams):
                if np.logical_and(beam >= 0, beam < max_beams):
                    beam_diff = beam_counter - 1
                    try:
                        current_beam_slist = fitacf_data[scan + beam_diff]['slist']
                    except (KeyError, IndexError):
                        continue
                    current_beam_gscat = fitacf_data[scan + beam_diff]['gflg']
                    isin = np.isin([gates[0], gates[1], gates[2]], current_beam_slist)
                    isin_indexes = np.where(isin)[0]
                    if isin_indexes.size != 0:
                        slist_indexes = np.where(np.isin(current_beam_slist, gates))[0]
                        gscat = np.where(current_beam_gscat[slist_indexes] == 0)
                        scatter[isin_indexes[gscat], beam_counter] = 1
                        vels.extend(fitacf_data[scan + beam_diff]['v'][slist_indexes[gscat]])
            cumulative_weight += np.sum(scatter * weighting_array[scan_counter])

    if cumulative_weight > weight_score:
        return np.median(vels)
    else:
        return []


def records(n_scans, n_beams=16, nrang=75, seed=0):
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n_scans):
        for beam in range(n_beams):
            slist = np.sort(rng.choice(nrang, rng.integers(0, nrang // 2), replace=False)).astype(np.int16)
            out.append({'bmnum': beam, 'nrang': nrang, 'slist': slist,
                        'gflg': (rng.random(slist.size) < 0.2).astype(np.int8),
                        'v': rng.normal(0, 500, slist.size).astype(np.float32)})
    return out


@pytest.mark.parametrize('seed', range(3))
def test_median_filter_record_matches_per_gate_filter(seed):
    data = records(5, seed=seed)

    for record in range(len(data)):
        medians, passed = median_filter_record(_WEIGHTING_ARRAY, data, record, 16)
        for gate in data[record]['slist']:
            expected = median_filter(_WEIGHTING_ARRAY, data, record, 16, gate)
            if isinstance(expected, list):
                assert not passed[gate]
            else:
                assert passed[gate]
                assert medians[gate] == expected
