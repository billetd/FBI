import cartopy.crs as ccrs
import numpy as np
from matplotlib import pyplot as plt, ticker, cm
from matplotlib.colors import Normalize
from FBI.grid import darn_grid_cell_index
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def plot_noon_line(apex, time, coord='mlt'):
    """
    Line from 80 to 90 degrees magnetic latitude, pointing towards magnetic noon
    :param apex: apexpy.Apex object
    :param time: datetime
    :param coord: Only 'mag' is supported
    """

    if coord == 'mag':
        mlon = (apex.mlt2mlon(12, time))
        plt.plot([mlon, mlon], [90, 80], transform=ccrs.PlateCarree(), color='m', zorder=1)

    else:
        print('Warning in plot_noon_line: coord can only be \'mag\'')


def plot_vecs_model_darn_grid(lompe, ax, coord='mag'):
    """
    Fit velocity vectors on the SuperDARN equal area grid, drawn thicker in cells with data
    :param lompe: dict - One record from readwrite.fbi_load_hdf5()
    :param ax: cartopy axis from axis.get_local_axis()
    :param coord: Only 'mag' draws anything
    :return: (thin quiver, thick quiver)
    """

    # Grid the data locations, so we can see where we want to highlight them when plotting
    data_cells = darn_grid_cell_index(lompe['mlons_los'], lompe['mlats_los'])

    # Get velocity vectors and points of grid
    v_emag = np.array(lompe['v_e_darngrid'])
    v_nmag = np.array(lompe['v_n_darngrid'])
    mlons = np.array(lompe['mlons_darngrid'])
    mlats = np.array(lompe['mlats_darngrid'])

    # Vectors in grid cells with data in them
    thick = np.isin(darn_grid_cell_index(mlons, mlats), data_cells[data_cells >= 0])

    # Vectors are plotted at the nearest whole degree of longitude, and not at the pole
    keep = ~(mlats > 89)
    mlons, mlats = np.round(mlons)[keep], mlats[keep]
    v_emag, v_nmag, thick = v_emag[keep], v_nmag[keep], thick[keep]

    # Rotate and scale vectors (https://github.com/SciTools/cartopy/issues/1179)
    u_src_crs = v_emag / np.cos(mlats / 180 * np.pi)
    v_src_crs = v_nmag
    magnitude = np.sqrt(v_emag ** 2 + v_nmag ** 2)
    magn_src_crs = np.sqrt(u_src_crs ** 2 + v_src_crs ** 2)

    # Put the stretched vectors back to the right length. Velocities are rounded to whole
    # m/s on write, so a cell can be zero in both components. Those get zero length.
    nonzero = magn_src_crs != 0
    u_plot = np.divide(u_src_crs * magnitude, magn_src_crs, out=np.zeros_like(magnitude), where=nonzero)
    v_plot = np.divide(v_src_crs * magnitude, magn_src_crs, out=np.zeros_like(magnitude), where=nonzero)

    colours_norm = Normalize(vmin=0, vmax=1000)
    if coord == 'mag':

        # Plot all the vectors at grid points normally
        quiv_thin = ax.quiver(mlons, mlats, u_plot, v_plot,
                              magnitude, norm=colours_norm, scale=2000, scale_units='inches', width=0.001,
                              headwidth=3, transform=ccrs.PlateCarree(), angles='xy', cmap='viridis', zorder=3)

        # Thicker vectors in the grid cells with data
        quiv_thick = ax.quiver(mlons[thick], mlats[thick], u_plot[thick], v_plot[thick],
                               magnitude[thick], norm=colours_norm, scale=2000, scale_units='inches', width=0.003,
                               headwidth=3, transform=ccrs.PlateCarree(), angles='xy', cmap='viridis', zorder=3)
    else:
        quiv_thick = None
        quiv_thin = None

    # Colour bar
    mappable = cm.ScalarMappable(norm=colours_norm, cmap='viridis')
    locator = ticker.MaxNLocator(symmetric=True, min_n_ticks=3, integer=True, nbins='auto')
    ticks = locator.tick_values(vmin=0, vmax=1000)

    # Add a small axis for the colorbar
    cax = inset_axes(ax,
                     width="90%",  # width of colorbar
                     height="3%",  # height of colorbar
                     loc='lower center',
                     bbox_to_anchor=(0, -0.05, 1, 1),
                     bbox_transform=ax.transAxes,
                     borderpad=0)
    cb = plt.colorbar(mappable, extend='max', ticks=ticks, cax=cax, orientation='horizontal')
    cb.set_label(r'Ionospheric Drift Velocity [ms$^{-1}$]')

    return quiv_thin, quiv_thick


def plot_potential_contours(lompe, ot, apex, time, coord='mag'):
    """
    Filled contours of the electric potential, with a colour bar
    :param lompe: dict - One record from readwrite.fbi_load_hdf5()
    :param ot: The projection of the axis for 'mag', or the polplot Polarplot for 'mlt'
    :param apex: apexpy.Apex object
    :param time: datetime, for magnetic local time
    :param coord: 'mag' or 'mlt'
    """

    V = np.array(lompe['e_pot_model'])/1000
    pot_mlat = np.array(lompe['mlats_model'])
    pot_mlon = np.array(lompe['mlons_model'])

    # Fixed range of contours, ten each side of zero
    vmax = 60
    pot_zmin = -vmax
    pot_zmax = vmax
    contour_spacing = int(np.floor(np.max([abs(pot_zmin), abs(pot_zmax)]) / 10))

    # Making the levels required, but skipping 0 as default to avoid a contour at 0 position (looks weird)
    contour_levels = [*range(pot_zmin, 0, contour_spacing), *range(contour_spacing,
                                                                   pot_zmax + contour_spacing, contour_spacing)]

    if coord == 'mag':
        # Convert to xy
        x, y, z = ot.transform_points(ccrs.PlateCarree(), pot_mlon, pot_mlat).T

        x_new = x[~np.isnan(x)]
        y_new = y[~np.isnan(x)]
        V_new = V[~np.isnan(x)]

        cs = plt.tricontourf(x_new, y_new, V_new.T, levels=contour_levels, zorder=2, cmap='RdBu',
                             vmax=pot_zmax, vmin=pot_zmin, extend='both', alpha=0.5)
    if coord == 'mlt':
        pot_mlt = (apex.mlon2mlt(np.array(pot_mlon), time)) * 15
        cs = ot.contourf(pot_mlat, pot_mlt / 15, V, levels=contour_levels, zorder=2, cmap='RdBu',
                         vmax=pot_zmax, vmin=pot_zmin, extend='both', alpha=0.5)

    # Colour bar
    locator = ticker.MaxNLocator(symmetric=True, min_n_ticks=3, integer=True, nbins='auto')
    ticks = locator.tick_values(vmin=pot_zmin, vmax=pot_zmax)
    cb = plt.colorbar(cs, extend='both', ticks=ticks)
    cb.set_label('Electric Potential [kV]')


def plot_data_locs(lompe, ax, apex=None, time=None, coord='mag'):
    """
    Locations of the SuperDARN data that went into the fit
    :param lompe: dict - One record from readwrite.fbi_load_hdf5()
    :param ax: cartopy axis for 'mag', or polplot Polarplot for 'mlt'
    :param apex: apexpy.Apex object, for 'mlt'
    :param time: datetime, for 'mlt'
    :param coord: 'mag' or 'mlt'
    """

    # Get coordinates
    data_mlats = lompe['mlats_los']
    data_mlons = lompe['mlons_los']

    if coord == 'mag':
        ax.scatter(data_mlons, data_mlats, s=0.4, color='k', zorder=1,
                   transform=ccrs.PlateCarree(), marker='x')
    elif coord == 'mlt':
        data_mlts = apex.mlon2mlt(data_mlons, time)
        ax.scatter(data_mlats, data_mlts, s=0.4, linewidth=0, color='black', zorder=1)


def plot_boundary_box(lompe, ax, apex, time):
    """
    Boundary of the fit, on a polar plot
    :param lompe: dict - One record from readwrite.fbi_load_hdf5()
    :param ax: polplot Polarplot
    :param apex: apexpy.Apex object
    :param time: datetime, for magnetic local time
    """

    # Get the coordinates of the boundary of the fit
    bound_mlts = apex.mlon2mlt(lompe['bound_mlons'], time)

    # Plot
    ax.plot(lompe['bound_mlats'], bound_mlts, color='black', linewidth=1, zorder=3)
