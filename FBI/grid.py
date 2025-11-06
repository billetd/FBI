"""
Code for handling grid related things, like making the lompe grid, or handling the equal area superdarn grid
"""
import lompe
import numpy as np
from geodarn.gridding import create_grid
from dataclasses import dataclass, field


@dataclass
class Container:
    location: np.ndarray = field(
        metadata={'group': 'data',
                  'units': 'degrees',
                  'description': 'Array of [lon, lat] locations'})


def sdarn_grid(apex):
    """

    :param apex:
    :param time:
    :return:
    """

    # Get the Superdarn grid
    _, _, darn_grid = create_grid(60, 1, 'north')
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
