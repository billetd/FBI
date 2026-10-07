import numpy as np
import pydarn
from FBI.los import Gates, los_azimuth


def test_azimuth_points_towards_the_radar():
    # Gates north, east, south and west of a radar on the equator
    azimuth = los_azimuth(0., 0., np.array([10., 0., -10., 0.]), np.array([0., 10., 0., -10.]))

    np.testing.assert_allclose(azimuth, [180, -90, 0, 90], atol=1e-5)


def test_gate_numbering():
    gates = Gates()
    inv = gates.table(pydarn.RadarID(64), 16, 75, 45, 180)
    rkn = gates.table(pydarn.RadarID(65), 16, 75, 45, 180)

    np.testing.assert_array_equal(inv, np.arange(1200).reshape(16, 75))
    np.testing.assert_array_equal(rkn, 1200 + np.arange(1200).reshape(16, 75))
    assert gates.table(pydarn.RadarID(64), 16, 75, 45, 180) is inv
    assert gates.size == 2400
    np.testing.assert_allclose(np.hypot(gates.le, gates.ln), 1)

    # Gates further down the beam are further from the radar
    radar = pydarn.SuperDARNRadars.radars[pydarn.RadarID(64)].hardware_info.geographic
    distance = np.hypot(gates.lat[inv[7]] - radar.lat, (gates.lon[inv[7]] - radar.lon) * np.cos(np.radians(radar.lat)))
    assert (np.diff(distance) > 0).all()
