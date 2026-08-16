"""
Functions for drawing individual plot elements (vectors, potential contours, data locations, ...)
onto an axis built by FBI.plotting.axis. Each function takes `record`, a single scan dict as
returned by FBI.readwrite.fbi_load_hdf5().

`ax` throughout accepts either a cartopy GeoAxes (for coord='mag', as returned by get_local_axis())
or a Polarplot instance (for coord='mlt', as returned by get_polar_axis()) - both expose compatible
plot/scatter/contourf-style methods.
"""
import cartopy.crs as ccrs
import numpy as np
from geodarn.gridding import create_grid_records
from matplotlib import ticker, cm
from matplotlib.colors import Normalize
from FBI.grid import Container
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

DEFAULT_VELOCITY_VMAX = 1000  # m/s
DEFAULT_POTENTIAL_VMAX_KV = 60  # kV


def _symmetric_ticks(vmin, vmax):
    locator = ticker.MaxNLocator(symmetric=True, min_n_ticks=3, integer=True, nbins='auto')
    return locator.tick_values(vmin=vmin, vmax=vmax)


def _underlying_axes(ax):
    """Resolve the real matplotlib Axes, whether `ax` is one already or a Polarplot wrapping one."""
    return getattr(ax, 'ax', ax)


def plot_noon_line(apex, time, ax, coord='mag', color='m', zorder=1):
    """
    Draw a line from the pole to the edge of the map indicating magnetic local noon.
    :param apex: apexpy.Apex instance.
    :param time: datetime of the scan.
    :param ax: cartopy GeoAxes to draw onto. Only coord='mag' is implemented.
    :param coord: 'mag' is the only supported value.
    """

    if coord != 'mag':
        raise NotImplementedError(f"plot_noon_line only supports coord='mag', got {coord!r}")

    mlon = apex.mlt2mlon(12, time)
    ax.plot([mlon, mlon], [90, 80], transform=ccrs.PlateCarree(), color=color, zorder=zorder)


def plot_vecs_model_darn_grid(record, ax, coord='mag', vmax=DEFAULT_VELOCITY_VMAX, cmap='viridis',
                              scale=2000, thin_width=0.001, thick_width=0.003, headwidth=3):
    """
    Plot the Lompe fit velocity vectors on the SuperDARN equal-area grid, thickening vectors at
    grid cells where line-of-sight data actually constrained the fit.
    :param record: scan dict, as returned by fbi_load_hdf5().
    :param ax: cartopy GeoAxes to draw onto. Only coord='mag' is implemented.
    :param coord: 'mag' is the only supported value.
    :param vmax: Upper bound of the velocity colour scale, in m/s.
    :return: (thin_quiver, thick_quiver)
    """

    if coord != 'mag':
        raise NotImplementedError(f"plot_vecs_model_darn_grid only supports coord='mag', got {coord!r}")

    # Get the locations of data
    data_mlats = record['mlats_los']
    data_mlons = record['mlons_los']

    # Grid the data locations, so we can see where we want to highlight them when plotting
    location = np.array([data_mlons, data_mlats]).T
    located = Container(location=location)
    idx_in_grid, darn_grid = create_grid_records(located)
    mlons_grid = darn_grid[:, 0]
    mlats_grid = darn_grid[:, 1]

    # Get velocity vectors and points of grid
    v_emag = np.array(record['v_e_darngrid'])
    v_nmag = np.array(record['v_n_darngrid'])
    mlons = np.array(record['mlons_darngrid'])

    # Seems to be a problem with floating point precision if this isn't done
    # This whole method could use re-working tbh. Currently quite janky.
    # TODO: Rework method for finding grid cells with data
    mlons = np.round(mlons)
    mlats = np.array(record['mlats_darngrid'])
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

    colours_norm = Normalize(vmin=0, vmax=vmax)

    # Plot all the vectors at grid points normally
    quiv_thin = ax.quiver(mlons, mlats, u_src_crs * magnitude / magn_src_crs, v_src_crs * magnitude / magn_src_crs,
                          magnitude, norm=colours_norm, scale=scale, scale_units='inches', width=thin_width,
                          headwidth=headwidth, transform=ccrs.PlateCarree(), angles='xy', cmap=cmap, zorder=3)

    # Figure out the velocities of data points which fall in the sdarn grid
    mlats_idx = mlats_grid[idx_in_grid].compressed()
    mlons_idx = mlons_grid[idx_in_grid].compressed()
    locs = []
    for this_mlat, this_mlon in zip(mlats_idx, mlons_idx):
        loc = np.where((this_mlat == mlats) & ((np.floor(this_mlon) == mlons) | (np.ceil(this_mlon) == mlons)))
        if len(loc[0]) > 0:
            locs.append(loc[0])

    u_src_crs_thick = u_src_crs[locs]
    v_src_crs_thick = v_src_crs[locs]
    magnitude_thick = magnitude[locs]
    magn_src_crs_thick = magn_src_crs[locs]
    thick_mlons = mlons[locs]
    thick_mlats = mlats[locs]

    # Plot thick vectors, to highlight grid cells actually constrained by data
    quiv_thick = ax.quiver(thick_mlons, thick_mlats, u_src_crs_thick * magnitude_thick /
                           magn_src_crs_thick, v_src_crs_thick * magnitude_thick / magn_src_crs_thick,
                           magnitude_thick, norm=colours_norm, scale=scale, scale_units='inches',
                           width=thick_width, headwidth=headwidth, transform=ccrs.PlateCarree(),
                           angles='xy', cmap=cmap, zorder=3)

    # Colour bar
    mappable = cm.ScalarMappable(norm=colours_norm, cmap=cmap)
    ticks = _symmetric_ticks(0, vmax)

    # Add a small axis for the colorbar
    cax = inset_axes(ax, width="90%", height="3%", loc='lower center', bbox_to_anchor=(0, -0.05, 1, 1),
                     bbox_transform=ax.transAxes, borderpad=0)
    cb = ax.figure.colorbar(mappable, extend='max', ticks=ticks, cax=cax, orientation='horizontal')
    cb.set_label(r'Ionospheric Drift Velocity [ms$^{-1}$]')

    return quiv_thin, quiv_thick


def plot_potential_contours(record, ax, apex, time, coord='mag', projection=None,
                            vmax_kv=DEFAULT_POTENTIAL_VMAX_KV, cmap='RdBu', alpha=0.5):
    """
    Contour-fill the Lompe fit electric potential.
    :param record: scan dict, as returned by fbi_load_hdf5().
    :param ax: cartopy GeoAxes for coord='mag', or a Polarplot for coord='mlt'.
    :param apex: apexpy.Apex instance, used to convert mlon to mlt for coord='mlt'.
    :param time: datetime of the scan.
    :param coord: 'mag' or 'mlt'.
    :param projection: Required for coord='mag'. The cartopy Projection (the `projection`/`ot`
        returned by get_local_axis()) used to convert mlat/mlon into the axis's projected coordinates.
    :param vmax_kv: The potential colour scale runs from -vmax_kv to +vmax_kv, in kV.
    :return: the contourf QuadContourSet
    """

    V = np.array(record['e_pot_model']) / 1000
    pot_mlat = np.array(record['mlats_model'])
    pot_mlon = np.array(record['mlons_model'])

    pot_zmin, pot_zmax = -vmax_kv, vmax_kv
    contour_spacing = max(int(np.floor(vmax_kv / 10)), 1)

    # Making the levels required, but skipping 0 as default to avoid a contour at 0 position (looks weird)
    contour_levels = [*range(pot_zmin, 0, contour_spacing),
                      *range(contour_spacing, pot_zmax + contour_spacing, contour_spacing)]

    if coord == 'mag':
        if projection is None:
            raise ValueError("plot_potential_contours requires `projection` when coord='mag'")

        # Convert to xy
        x, y, _ = projection.transform_points(ccrs.PlateCarree(), pot_mlon, pot_mlat).T
        finite = ~np.isnan(x)
        cs = ax.tricontourf(x[finite], y[finite], V[finite], levels=contour_levels, zorder=2, cmap=cmap,
                            vmax=pot_zmax, vmin=pot_zmin, extend='both', alpha=alpha)
    elif coord == 'mlt':
        pot_mlt = apex.mlon2mlt(pot_mlon, time) * 15
        cs = ax.contourf(pot_mlat, pot_mlt / 15, V, levels=contour_levels, zorder=2, cmap=cmap,
                         vmax=pot_zmax, vmin=pot_zmin, extend='both', alpha=alpha)
    else:
        raise NotImplementedError(f"plot_potential_contours only supports coord in ('mag', 'mlt'), got {coord!r}")

    # Colour bar
    real_ax = _underlying_axes(ax)
    ticks = _symmetric_ticks(pot_zmin, pot_zmax)
    cb = real_ax.figure.colorbar(cs, ax=real_ax, extend='both', ticks=ticks)
    cb.set_label('Electric Potential [kV]')

    return cs


def plot_data_locs(record, ax, apex=None, time=None, coord='mag', size=0.4, color='k'):
    """
    Scatter the locations of the line-of-sight data going into the fit.
    :param ax: cartopy GeoAxes for coord='mag', or a Polarplot for coord='mlt'.
    :param apex: apexpy.Apex instance. Required for coord='mlt'.
    :param time: datetime of the scan. Required for coord='mlt'.
    """

    data_mlats = record['mlats_los']
    data_mlons = record['mlons_los']

    if coord == 'mag':
        ax.scatter(data_mlons, data_mlats, s=size, color=color, zorder=1,
                  transform=ccrs.PlateCarree(), marker='x')
    elif coord == 'mlt':
        data_mlts = apex.mlon2mlt(data_mlons, time)
        ax.scatter(data_mlats, data_mlts, s=size, linewidth=0, color=color, zorder=1)
    else:
        raise NotImplementedError(f"plot_data_locs only supports coord in ('mag', 'mlt'), got {coord!r}")


def plot_boundary_box(record, ax, apex, time, color='black', linewidth=1):
    """
    Draw the outer boundary of the Lompe fit grid. Only meaningful on a Polarplot (coord='mlt') axis.
    """

    bound_mlts = apex.mlon2mlt(record['bound_mlons'], time)
    ax.plot(record['bound_mlats'], bound_mlts, color=color, linewidth=linewidth, zorder=3)
