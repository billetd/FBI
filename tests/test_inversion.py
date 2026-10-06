"""
FBI.inversion is a cut down copy of lompe's Emodel and run_inversion(), so check that it still
gives the same matrices and model vector as stock lompe.
"""
import numpy as np
import pytest
import lompe
from FBI.inversion import Model


@pytest.fixture(scope='module')
def grid():
    projection = lompe.cs.CSprojection((-100, 62), -30)
    return lompe.cs.CSgrid(projection, 1500e3, 1200e3, 100e3, 100e3, R=6481.2e3)


@pytest.fixture(scope='module')
def emodel(grid):
    return lompe.Emodel(grid, Hall_Pedersen_conductance=None, ew_regularization_limit=(50, 75))


@pytest.fixture(scope='module')
def model(grid):
    return Model(grid, l1=10, l2=0.1, ew_regularization_limit=(50, 75))


def los_data(n, seed=0):
    """Random line-of-sight data around the grid, with a few points outside the data density grid"""
    rng = np.random.default_rng(seed)
    lon = rng.uniform(-125, -75, n)
    lat = rng.uniform(52, 72, n)
    lat[:5] = rng.uniform(20, 25, 5)  # Well outside the grid
    vlos = np.abs(rng.normal(0, 500, n))
    azimuth = rng.uniform(0, 2 * np.pi, n)
    le, ln = np.sin(azimuth), np.cos(azimuth)
    error = rng.uniform(20, 100, n)
    return lon, lat, vlos, le, ln, error


def test_same_matrices_as_emodel(model, emodel):
    np.testing.assert_array_equal(model.lat_E, emodel.lat_E)
    np.testing.assert_array_equal(model.lon_E, emodel.lon_E)
    for ours, lompes in zip(model.v_matrix(), (emodel.Ve, emodel.Vn)):
        np.testing.assert_array_equal(ours, lompes)

    lon, lat = los_data(50)[:2]
    for ours, lompes in zip(model.v_matrix(lon, lat), emodel._v_matrix(lon=lon, lat=lat)):
        np.testing.assert_array_equal(ours, lompes)

    LTL = 10 * np.eye(emodel.grid_E.size) + 0.1 * emodel.LTLe / np.median(emodel.LTLe.diagonal())
    np.testing.assert_array_equal(model.LTL, LTL)


def test_matches_stock_run_inversion(model, emodel):
    lon, lat, vlos, le, ln, error = los_data(400)

    used = model.biggrid.ingrid(lon, lat)
    m = model.solve(model.los_matrix(lon[used], lat[used], le[used], ln[used]),
                    lon[used], lat[used], vlos[used], error[used])

    emodel.clear_model()
    emodel.add_data(lompe.Data(vlos, coordinates=np.vstack((lon, lat)), LOS=np.vstack((le, ln)),
                               datatype='convection', error=error, iweight=1.0))
    emodel.run_inversion(l1=10, l2=0.1, lapack_driver='gelsy')

    # Same data points used
    assert not used[:5].any()
    np.testing.assert_array_equal(emodel.data['convection'][0].coords['lon'], lon[used])
    np.testing.assert_array_equal(emodel.data['convection'][0].coords['lat'], lat[used])

    # Same model vector
    np.testing.assert_allclose(m, emodel.m, rtol=1e-10, atol=1e-10 * np.max(np.abs(emodel.m)))
