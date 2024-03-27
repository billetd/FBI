import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
import numpy as np
import matplotlib.path as mpltPath
from matplotlib import pyplot as plt
from pydarn.plotting.projections import convert_geo_coastline_to_mag
from shapely.geometry import MultiLineString


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
    gl.xlines = False

    # Read in the geometry object of the coastlines
    cc = cfeature.NaturalEarthFeature('physical', 'coastline', '50m',
                                      color='k', zorder=2.0)

    # Plot coastlines
    for shape in list(cc.geometries()):
        if isinstance(shape, MultiLineString):  # Don't plot multi geoms as it breaks
            continue
        glats = shape.coords.xy[1]
        glons = shape.coords.xy[0]
        mlats, mlons = apex.geo2apex(glats, glons, 300)
        x_coast, y_coast, zcoast = ot.transform_points(ccrs.PlateCarree(), mlons, mlats).T
        points = path.contains_points(list(zip(x_coast, y_coast)))
        if any(points):  # Check if any of the points to plot are actually in the plot window
            # plt.fill(x_coast, y_coast, zorder=0, color='grey')  # Doesn't work right atm. Weird shapes.
            plt.plot(x_coast, y_coast, zorder=0, color='grey', linewidth=0.5, alpha=0.6)

    return ax, ot, 'mag', fig
