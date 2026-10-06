"""
Code for handling grid related things, like making the lompe grid, or handling the equal area superdarn grid
"""
import lompe
import numpy as np
from geodarn.gridding import create_grid

# Lower latitude boundary and cell height of the SuperDARN equal area grid [degrees]
DARN_GRID_LAT_MIN = 60
DARN_GRID_LAT_WIDTH = 1


def sdarn_grid(apex):
    """

    :param apex:
    :param time:
    :return:
    """

    # Get the Superdarn grid
    _, _, darn_grid = create_grid(DARN_GRID_LAT_MIN, DARN_GRID_LAT_WIDTH, 'north')
    # _, _, darn_grid = create_grid(10, 1, 'north')
    darn_grid = darn_grid.reshape(-1, 2)
    mlats_darngrid = darn_grid[:, 1].compressed()
    mlons_darngrid = darn_grid[:, 0].compressed()

    glats_darngrid, glons_darngrid, _ = apex.apex2geo(mlats_darngrid, mlons_darngrid, 300)
    darn_grid_stuff = {'mlats_darngrid': mlats_darngrid,
                       'mlons_darngrid': mlons_darngrid,
                       'glats_darngrid': glats_darngrid,
                       'glons_darngrid': glons_darngrid}

    return darn_grid_stuff


def darn_grid_cell_index(mlons, mlats, lat_min=DARN_GRID_LAT_MIN, lat_width=DARN_GRID_LAT_WIDTH):
    """
    Find which cell of the SuperDARN equal area grid each point is in
    :param mlons: array - Magnetic longitudes [degrees]
    :param mlats: array - Magnetic latitudes [degrees]
    :param lat_min: Lower latitude boundary of the grid [degrees]
    :param lat_width: Cell height [degrees]
    :return: array of int - Index of each point's cell in the flattened geodarn
             create_grid() grid, or -1 for points outside the grid
    """

    lat_divs_bottom, lon_divs, darn_grid = create_grid(lat_min, lat_width, 'north')
    max_num_lons = darn_grid.shape[1]

    mlons = (np.asarray(mlons, dtype=float) + 180) % 360 - 180
    mlats = np.asarray(mlats, dtype=float)
    cells = np.full(mlats.shape, -1)

    # Latitude band of each point, then the longitude cell within that band
    on_grid = np.isfinite(mlons) & (mlats >= lat_min) & (mlats <= 90)
    rows = np.full(mlats.shape, -1)
    rows[on_grid] = np.minimum(np.searchsorted(lat_divs_bottom, mlats[on_grid], side='right') - 1,
                               len(lat_divs_bottom) - 1)
    for row in np.unique(rows[on_grid]):
        in_row = rows == row
        cols = np.searchsorted(lon_divs[row], mlons[in_row], side='right') - 1
        cells[in_row] = row * max_num_lons + np.clip(cols, 0, len(lon_divs[row]) - 2)

    return cells


def lompe_grid_canada(apex):
    """
    Creates a lompe grid ideal for doing lompe maps with SuperDARN Canada Borealis radars
    :param apex:
    :return:
    """

    # cubed sphere grid parameters:
    mag_position = (-34, 81)
    # mag_position = (-34, 70)
    lat, lon, z = apex.apex2geo(mag_position[1], mag_position[0], 300)
    position = (lon, lat)  # lon, lat for center of the grid

    # Current one
    orientation = -127.7
    l, w, lres, wres = 5000e3, 4000e3, 75.e3, 75.e3  # slightly shorter, for polar plotting
    # l, w, lres, wres = 16000e3, 16000e3, 150.e3, 150.e3  # for bill
    # l, w, lres, wres = 10000e3, 10000e3, 150.e3, 150.e3

    # Create grid object:
    grid = lompe.cs.CSgrid(lompe.cs.CSprojection(position, orientation), l, w, lres, wres, R=6481.2e3)

    return grid
