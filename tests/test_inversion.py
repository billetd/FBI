"""
FBI.inversion is a cut down copy of lompe.Emodel.run_inversion(), so check that it still
gives the same model vector as stock lompe.
"""
import numpy as np
import pytest
import lompe
from FBI.inversion import prepare_inversion, solve_los


@pytest.fixture(scope='module')
def model():
    projection = lompe.cs.CSprojection((-100, 62), -30)
    grid = lompe.cs.CSgrid(projection, 1500e3, 1200e3, 100e3, 100e3, R=6481.2e3)
    return lompe.Emodel(grid, Hall_Pedersen_conductance=None, ew_regularization_limit=(50, 75))


def los_data(n, seed=0):
    """Random line-of-sight data around the grid, with a few points outside the data density grid"""
    rng = np.random.default_rng(seed)
    lon = rng.uniform(-125, -75, n)
    lat = rng.uniform(52, 72, n)
    lat[:5] = rng.uniform(20, 25, 5)  # Well outside the grid
    vlos = np.abs(rng.normal(0, 500, n))
    vlos[5] = np.nan
    azimuth = rng.uniform(0, 2 * np.pi, n)
    le, ln = np.sin(azimuth), np.cos(azimuth)
    error = rng.uniform(20, 100, n)
    return lon, lat, vlos, le, ln, error


def test_matches_stock_run_inversion(model):
    lon, lat, vlos, le, ln, error = los_data(400)

    m, used = solve_los(model, prepare_inversion(model), lon, lat, vlos, le, ln, error)

    model.clear_model()
    model.add_data(lompe.Data(vlos, coordinates=np.vstack((lon, lat)), LOS=np.vstack((le, ln)),
                              datatype='convection', error=error, iweight=1.0))
    model.run_inversion(l1=10, l2=0.1, lapack_driver='gelsy')

    # Same data points used
    assert not used[:6].any()
    np.testing.assert_array_equal(model.data['convection'][0].coords['lon'], lon[used])
    np.testing.assert_array_equal(model.data['convection'][0].coords['lat'], lat[used])

    # Same model vector
    np.testing.assert_allclose(m, model.m, rtol=1e-10, atol=1e-10 * np.max(np.abs(model.m)))


def test_too_little_data(model):
    lon, lat, vlos, le, ln, error = los_data(10)

    # Only the points outside the grid, or a single point inside it
    for keep in (slice(0, 5), slice(0, 7)):
        m, used = solve_los(model, prepare_inversion(model), lon[keep], lat[keep], vlos[keep],
                            le[keep], ln[keep], error[keep])
        assert m is None
        assert used.sum() <= 1
