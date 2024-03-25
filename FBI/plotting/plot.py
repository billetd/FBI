import cartopy.crs as ccrs
import numpy as np
from geodarn.gridding import create_grid_records
from matplotlib import pyplot as plt, ticker, cm
from matplotlib.colors import Normalize
from FBI.grid import Container


def plot_noon_line(apex, time, coord='mlt'):
    """

    :param apex:
    :param time:
    :param coord:
    :return:
    """

    if coord == 'mag':
        mlon = (apex.mlt2mlon(12, time))
        plt.plot([mlon, mlon], [90, 80], transform=ccrs.PlateCarree(), color='m', zorder=1)

    else:
        print('Warning in plot_noon_line: coord can only be \'mag\'')


def plot_vecs_model_darn_grid(lompe, ax, mlon=None):
    """

    :param lompe:
    :param ax:
    :param mlon:
    :return:
    """

    # Get the locations of data
    data_mlats = lompe['mlats_los']
    data_mlons = lompe['mlons_los']

    # Grid the data locations, so we can see where we want to highlight them when plotting
    # Put in a dictionary for the gridding code
    # from BorealisConvection.lompe_gridding import Container
    location = np.array([data_mlons, data_mlats]).T
    located = Container(location=location)
    idx_in_grid, darn_grid = create_grid_records(located)
    mlons_grid = darn_grid[:, 0]
    mlats_grid = darn_grid[:, 1]

    # Get velocity vectors and points of grid
    v_emag = np.array(lompe['v_e_darngrid'])
    v_nmag = np.array(lompe['v_n_darngrid'])
    mlons = np.array(lompe['mlons_darngrid'])

    # Seems to be a problem with floating point precision if this isn't done
    # This whole method could use re-working tbh. Currently quite janky.
    # TODO: Rework method for finding grid cells with data
    mlons = np.round(mlons)
    mlats = np.array(lompe['mlats_darngrid'])
    hilats = np.where(mlats > 89)
    mlats = np.delete(mlats, hilats)
    mlons = np.delete(mlons, hilats)
    v_emag = np.delete(v_emag, hilats)
    v_nmag = np.delete(v_nmag, hilats)

    # Rotate and scale vectors (https://github.com/SciTools/cartopy/issues/1179)
    u_src_crs = v_emag / np.cos(mlats / 180 * np.pi)
    v_src_crs = v_nmag
    magnitude = np.sqrt(v_emag ** 2 + v_nmag ** 2)
    magn_src_crs = np.sqrt(u_src_crs ** 2 + v_src_crs ** 2)

    colours_norm = Normalize(vmin=0, vmax=1000)
    if mlon is True:

        # Plot all the vectors at grid points normally
        quiv_thin = ax.quiver(mlons, mlats, u_src_crs * magnitude / magn_src_crs, v_src_crs * magnitude / magn_src_crs,
                              magnitude, norm=colours_norm, scale=2000, scale_units='inches', width=0.001,
                              headwidth=3, transform=ccrs.PlateCarree(), angles='xy', cmap='viridis', zorder=3)

        # Figure out the velocities of data points which fall in the sdarn grid
        mlats_idx = mlats_grid[idx_in_grid].compressed()
        mlons_idx = mlons_grid[idx_in_grid].compressed()
        locs = []
        for this_mlat, this_mlon in zip(mlats_idx, mlons_idx):
            loc = np.where((this_mlat == mlats) & ((np.floor(this_mlon) == mlons) | (np.ceil(this_mlon) == mlons)))
            if loc[0]:
                locs.append(loc[0])

        u_src_crs_thick = u_src_crs[locs]
        v_src_crs_thick = v_src_crs[locs]
        magnitude_thick = magnitude[locs]
        magn_src_crs_thick = magn_src_crs[locs]
        thick_mlons = mlons[locs]
        thick_mlats = mlats[locs]

        # Plot thick vectors
        quiv_thick = ax.quiver(thick_mlons, thick_mlats, u_src_crs_thick * magnitude_thick /
                               magn_src_crs_thick, v_src_crs_thick * magnitude_thick / magn_src_crs_thick,
                               magnitude_thick, norm=colours_norm, scale=2000, scale_units='inches', width=0.003,
                               headwidth=3, transform=ccrs.PlateCarree(), angles='xy', cmap='viridis', zorder=3)
    else:
        quiv_thick = None
        quiv_thin = None

    # Colour bar
    mappable = cm.ScalarMappable(norm=colours_norm, cmap='viridis')
    locator = ticker.MaxNLocator(symmetric=True, min_n_ticks=3, integer=True, nbins='auto')
    ticks = locator.tick_values(vmin=0, vmax=1000)

    # Add a small axis for the colorbar
    sub_ax = plt.axes([0.77, 0.1, 0.03, 0.8])
    cb = plt.colorbar(mappable, extend='max', ticks=ticks, cax=sub_ax)
    cb.set_label(r'Ionospheric Drift Velocity [ms$^{-1}$]')

    return quiv_thin, quiv_thick
