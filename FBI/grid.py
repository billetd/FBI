"""
Code for handling grid related things, like making the lompe grid, or handling the equal area superdarn grid
"""
import numpy as np
from secsy import cubedsphere as cs

# Lower latitude boundary and cell height of the SuperDARN equal area grid [degrees]
DARN_GRID_LAT_MIN = 60
DARN_GRID_LAT_WIDTH = 1


def equal_area_grid(lat_min=DARN_GRID_LAT_MIN, lat_width=DARN_GRID_LAT_WIDTH):
    """
    The SuperDARN equal area grid in the northern hemisphere, the same as geodarn's create_grid()
    :param lat_min: Lower latitude boundary of the grid [degrees]
    :param lat_width: Cell height [degrees]
    :return: (lat_divs_bottom, lon_divs) - The lower latitude of each row of cells, and the
             longitude boundaries of the cells in each row
    """

    lat_divs_bottom = np.arange(lat_min, 90, lat_width)
    cos_lat = np.cos(np.deg2rad(lat_divs_bottom + lat_width / 2))
    num_lons = np.int32(np.rint(360.0 / (lat_width / cos_lat)))
    lon_divs = [np.linspace(-180, 180, n + 1) for n in num_lons]

    return lat_divs_bottom, lon_divs


def sdarn_grid(apex):
    """
    Centres of the SuperDARN equal area grid cells
    :param apex: apexpy.Apex object
    :return: dict of the magnetic and geographic coordinates of the cells, in the order
             darn_grid_cell_index() numbers them
    """

    lat_divs_bottom, lon_divs = equal_area_grid()
    lat_centers = lat_divs_bottom + DARN_GRID_LAT_WIDTH / 2
    mlons_darngrid = np.concatenate([divs[:-1] + (divs[1] - divs[0]) / 2 for divs in lon_divs])
    mlats_darngrid = np.concatenate([np.full(divs.size - 1, lat) for divs, lat in zip(lon_divs, lat_centers)])

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
    :return: array of int - Index of each point's cell, counting along each row from the lowest,
             or -1 for points outside the grid
    """

    lat_divs_bottom, lon_divs = equal_area_grid(lat_min, lat_width)
    row_start = np.cumsum([0] + [divs.size - 1 for divs in lon_divs])

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
        cells[in_row] = row_start[row] + np.clip(cols, 0, len(lon_divs[row]) - 2)

    return cells


def lompe_grid_canada(apex):
    """
    Creates a lompe grid ideal for doing lompe maps with SuperDARN Canada Borealis radars
    :param apex: apexpy.Apex object
    :return: secsy CSgrid
    """

    # Centre of the grid, from its magnetic position
    lat, lon, _ = apex.apex2geo(81, -34, 300)

    # Slightly shorter, for polar plotting
    l, w, lres, wres = 5000e3, 4000e3, 75.e3, 75.e3

    return cs.CSgrid(cs.CSprojection((lon, lat), -127.7), l, w, lres, wres, R=6481.2e3)
