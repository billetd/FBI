import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
import numpy as np
import matplotlib.path as mpltPath
from matplotlib import pyplot as plt
from polplot import Polarplot
from shapely.geometry import MultiLineString

# Coastlines in projected coordinates, keyed on the apex epoch, projection and window
_coastline_cache = {}


def coastline_segments(apex, ot, path, window):
    """
    Coastline segments in magnetic coordinates, projected and ready to plot
    :param apex: apexpy.Apex object
    :param ot: cartopy projection of the axis
    :param path: matplotlib Path of the plot window
    :param window: plot window as (x0, x1, y0, y1), for the cache key
    :return: list of (x, y) arrays, one per segment in view
    """

    key = (apex.year, ot.proj4_init, window)
    segments = _coastline_cache.get(key)
    if segments is not None:
        return segments

    # Read in the geometry object of the coastlines
    cc = cfeature.NaturalEarthFeature('physical', 'coastline', '50m',
                                      color='k', zorder=2.0)
    shapes = [shape.coords.xy for shape in cc.geometries()
              if not isinstance(shape, MultiLineString)]  # Don't plot multi geoms as it breaks

    # All the points in one apex call and one projection, rather than one per shape
    glons = np.concatenate([coords[0] for coords in shapes])
    glats = np.concatenate([coords[1] for coords in shapes])
    ends = np.cumsum([len(coords[0]) for coords in shapes])

    mlats, mlons = apex.geo2apex(glats, glons, 300)
    x_coast, y_coast, zcoast = ot.transform_points(ccrs.PlateCarree(), mlons, mlats).T

    # Keep segments with any point in the plot window
    inside = path.contains_points(np.column_stack((x_coast, y_coast)))
    segments = [(x_coast[start:end], y_coast[start:end])
                for start, end in zip(np.concatenate(([0], ends[:-1])), ends)
                if inside[start:end].any()]

    _coastline_cache[key] = segments
    return segments


def get_local_axis(apex):
    """

    :return:
    """

    # Set the projection to orthographic
    ot = ccrs.Orthographic(-61.15, 83)
    fig, ax = plt.subplots(subplot_kw={'projection': ot}, gridspec_kw={'wspace': 0.05, 'hspace': 0.05})

    # Set up the plot
    pos_lower = [-82, 56]
    pos_higher = [64, 71]
    # pos_lower = [-82, 40]
    # pos_higher = [50, 50]
    xs, ys, zs = ot.transform_points(ccrs.PlateCarree(), np.array((pos_lower[0], pos_lower[1])),
                                     np.array((pos_higher[0], pos_higher[1]))).T

    # For checking if plots are in the plot window later
    path = mpltPath.Path([[xs[0], ys[0]], [xs[0], ys[1]], [xs[1], ys[1]], [xs[1], ys[0]]])

    # Setup the axis to look nice
    ax.set_xlim(xs)
    ax.set_ylim(ys)
    gl = ax.gridlines(draw_labels=True, y_inline=True, crs=ccrs.PlateCarree(), zorder=5, color='grey')
    gl.ylocator = mticker.FixedLocator(np.arange(50, 91, 10))
    gl.ylabel_style = {'color': 'black'}
    gl.xlocator = mticker.FixedLocator(np.arange(-180, 181, 45))
    gl.top_labels = False
    gl.bottom_labels = False
    gl.geo_labels = False
    gl.xlines = False

    # Plot coastlines
    for x_coast, y_coast in coastline_segments(apex, ot, path, (xs[0], xs[1], ys[0], ys[1])):
        # plt.fill(x_coast, y_coast, zorder=0, color='grey')  # Doesn't work right atm. Weird shapes.
        plt.plot(x_coast, y_coast, zorder=0, color='grey', linewidth=0.5, alpha=0.6)
    return ax, ot, 'mag', fig


def get_polar_axis(time, apex):

    fig = plt.figure()
    rect = [0.1, 0.1, 0.8, 0.8]
    ax = fig.add_axes(rect)

    pax = Polarplot(ax, minlat=50, linewidth=0.7)
    pax.coastlines(time=time, mag=apex, linewidth=0.5, resolution='50m', color='grey', alpha=0.5)

    lowlat_mlts = np.linspace(0, 24, num=360)
    lowlat_lats = np.zeros(360) + 50

    # I don't like the dashed line at the low lat boundary, so I'll just draw a solid one
    pax.plot(lowlat_lats, lowlat_mlts, color='black')

    return pax, 'mlt', fig
