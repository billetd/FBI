import datetime as dt
import numpy as np
import apexpy
import pytest
from FBI.grid import lompe_grid_canada
from FBI.readwrite import _to_qd


@pytest.fixture(scope='module')
def setup():
    apex = apexpy.Apex(dt.datetime(2025, 2, 24, 18), refh=300)
    grid = lompe_grid_canada(apex)
    rng = np.random.default_rng(0)
    i = rng.choice(grid.size, 300, replace=False)
    glat, glon = grid.lat.ravel()[i], grid.lon.ravel()[i]
    azimuth = rng.uniform(0, 2 * np.pi, glat.size)
    speed = rng.uniform(100, 2000, glat.size)
    v_e, v_n = speed * np.sin(azimuth), speed * np.cos(azimuth)
    return apex, glat, glon, v_e, v_n


def test_direction_matches_qd_grid(setup):
    apex, glat, glon, v_e, v_n = setup
    v_e_mag, v_n_mag = _to_qd(apex.basevectors_qd(glat, glon, 300, coords='geo'), v_e, v_n)

    # Step a short way along each vector and see where it goes in QD coordinates
    r = 6371.2 + 300
    step = 0.5 / np.hypot(v_e, v_n)
    qlat0, qlon0 = apex.geo2qd(glat, glon, 300)
    qlat1, qlon1 = apex.geo2qd(glat + np.degrees(step * v_n / r),
                               glon + np.degrees(step * v_e / (r * np.cos(np.radians(glat)))), 300)
    dq_e = ((qlon1 - qlon0 + 180) % 360 - 180) * np.cos(np.radians(qlat0))
    dq_n = qlat1 - qlat0

    angle = np.degrees(np.arctan2(dq_e * v_n_mag - dq_n * v_e_mag, dq_e * v_e_mag + dq_n * v_n_mag))
    assert np.abs(angle).max() < 1


def test_speed_is_kept(setup):
    apex, glat, glon, v_e, v_n = setup
    v_e_mag, v_n_mag = _to_qd(apex.basevectors_qd(glat, glon, 300, coords='geo'), v_e, v_n)

    np.testing.assert_allclose(np.hypot(v_e_mag, v_n_mag), np.hypot(v_e, v_n), rtol=1e-12)
    np.testing.assert_array_equal(_to_qd(apex.basevectors_qd(glat, glon, 300, coords='geo'),
                                         np.zeros(glat.size), np.zeros(glat.size)), 0)


def test_along_qd_meridian(setup):
    apex, glat, glon, _, _ = setup
    f1, f2 = apex.basevectors_qd(glat, glon, 300, coords='geo')

    # Moving along f2 doesn't change QD longitude
    v_e_mag, v_n_mag = _to_qd((f1, f2), 1000 * f2[0], 1000 * f2[1])

    np.testing.assert_allclose(v_e_mag, 0, atol=1e-9)
    assert (v_n_mag > 0).all()
