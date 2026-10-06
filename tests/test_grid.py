import datetime as dt
import apexpy
import numpy as np
import pytest
from FBI.grid import darn_grid_cell_index, equal_area_grid, sdarn_grid, DARN_GRID_LAT_MIN, DARN_GRID_LAT_WIDTH

lat_divs_bottom, lon_divs = equal_area_grid()
row_start = np.cumsum([0] + [divs.size - 1 for divs in lon_divs])


def test_same_grid_as_geodarn():
    gridding = pytest.importorskip('geodarn.gridding')
    lat_divs_geodarn, lon_divs_geodarn, _ = gridding.create_grid(DARN_GRID_LAT_MIN, DARN_GRID_LAT_WIDTH, 'north')

    np.testing.assert_array_equal(lat_divs_bottom, lat_divs_geodarn)
    assert len(lon_divs) == len(lon_divs_geodarn)
    for ours, geodarns in zip(lon_divs, lon_divs_geodarn):
        np.testing.assert_array_equal(ours, geodarns)


def test_points_are_inside_their_cell():
    rng = np.random.default_rng(0)
    mlons = rng.uniform(-180, 180, 20000)
    mlats = rng.uniform(DARN_GRID_LAT_MIN, 90, 20000)

    cells = darn_grid_cell_index(mlons, mlats)

    assert (cells >= 0).all()
    rows = np.searchsorted(row_start, cells, side='right') - 1
    cols = cells - row_start[rows]
    assert (lat_divs_bottom[rows] <= mlats).all()
    assert (mlats < lat_divs_bottom[rows] + DARN_GRID_LAT_WIDTH).all()
    lon_lo = np.array([lon_divs[r][c] for r, c in zip(rows, cols)])
    lon_hi = np.array([lon_divs[r][c + 1] for r, c in zip(rows, cols)])
    assert (lon_lo <= mlons).all() and (mlons < lon_hi).all()


def test_cell_centres_map_to_their_own_cell():
    darn_grid = sdarn_grid(apexpy.Apex(dt.datetime(2025, 2, 24), refh=300))

    cells = darn_grid_cell_index(darn_grid['mlons_darngrid'], darn_grid['mlats_darngrid'])

    np.testing.assert_array_equal(cells, np.arange(row_start[-1]))


def test_points_in_a_single_cell():
    # geodarn's create_grid_records() misses these, as no cell boundary lies between them
    cells = darn_grid_cell_index(np.array([-100.3, -100.2]), np.array([75.5, 75.6]))

    assert (cells >= 0).all() and cells[0] == cells[1]


def test_points_off_the_grid():
    cells = darn_grid_cell_index(np.array([10., np.nan, 10.]), np.array([DARN_GRID_LAT_MIN - 0.5, 70., np.nan]))

    np.testing.assert_array_equal(cells, -1)


def test_longitude_wraps():
    cells = darn_grid_cell_index(np.array([180., -180., 540.]), np.array([70.5, 70.5, 70.5]))

    assert cells[0] == cells[1] == cells[2] == row_start[10]
