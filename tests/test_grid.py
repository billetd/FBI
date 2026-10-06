import numpy as np
from geodarn.gridding import create_grid
from FBI.grid import darn_grid_cell_index, DARN_GRID_LAT_MIN, DARN_GRID_LAT_WIDTH

lat_divs_bottom, lon_divs, darn_grid = create_grid(DARN_GRID_LAT_MIN, DARN_GRID_LAT_WIDTH, 'north')
max_num_lons = darn_grid.shape[1]


def test_points_are_inside_their_cell():
    rng = np.random.default_rng(0)
    mlons = rng.uniform(-180, 180, 20000)
    mlats = rng.uniform(DARN_GRID_LAT_MIN, 90, 20000)

    cells = darn_grid_cell_index(mlons, mlats)

    assert (cells >= 0).all()
    rows, cols = cells // max_num_lons, cells % max_num_lons
    assert (lat_divs_bottom[rows] <= mlats).all()
    assert (mlats < lat_divs_bottom[rows] + DARN_GRID_LAT_WIDTH).all()
    lon_lo = np.array([lon_divs[r][c] for r, c in zip(rows, cols)])
    lon_hi = np.array([lon_divs[r][c + 1] for r, c in zip(rows, cols)])
    assert (lon_lo <= mlons).all() and (mlons < lon_hi).all()
    assert not darn_grid.reshape(-1, 2).mask[cells].any()


def test_cell_centres_map_to_their_own_cell():
    flat = darn_grid.reshape(-1, 2)
    # create_grid() leaves one unmasked (0, 0) entry at the end of each row, which isn't a cell
    on_grid = ~flat.mask[:, 0] & (flat.data[:, 1] >= DARN_GRID_LAT_MIN)
    centres = flat[on_grid].data

    cells = darn_grid_cell_index(centres[:, 0], centres[:, 1])

    np.testing.assert_array_equal(cells, np.flatnonzero(on_grid))


def test_points_in_a_single_cell():
    # geodarn's create_grid_records() misses these, as no cell boundary lies between them
    cells = darn_grid_cell_index(np.array([-100.3, -100.2]), np.array([75.5, 75.6]))

    assert (cells >= 0).all() and cells[0] == cells[1]


def test_points_off_the_grid():
    cells = darn_grid_cell_index(np.array([10., np.nan, 10.]), np.array([DARN_GRID_LAT_MIN - 0.5, 70., np.nan]))

    np.testing.assert_array_equal(cells, -1)


def test_longitude_wraps():
    cells = darn_grid_cell_index(np.array([180., -180., 540.]), np.array([70.5, 70.5, 70.5]))

    assert cells[0] == cells[1] == cells[2] == 10 * max_num_lons
