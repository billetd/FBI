import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
from pydarn.plotting.projections import convert_geo_coastline_to_mag


def get_local_axis(time):
    """

    :param time:
    :return:
    """

    # Set the projection to orthographic
    ot = ccrs.Orthographic(-61.15, 83)
    fig, ax = plt.subplots(subplot_kw={'projection': ot}, gridspec_kw={'wspace': 0.05, 'hspace': 0.05})

    # Set up the plot
    xs, ys, zs = ot.transform_points(ccrs.PlateCarree(), np.array((-82, 56)), np.array((64, 71))).T

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

    # Convert geometry object coordinates to mlat-mlt
    geom_mag = []
    for geom in cc.geometries():
        thing = convert_geo_coastline_to_mag(geom, time, mag_lon=True)
        if thing is not None:
            geom_mag.append(thing)
    cc_mag = cfeature.ShapelyFeature(geom_mag, ccrs.PlateCarree(), color='k', zorder=2.0)

    # Plot each geometry object
    for geom in cc_mag.geometries():
        mlon, mlat = np.array(geom.coords.xy[0]), np.array(geom.coords.xy[1])

        # Convert to x and y for plotting because it can't handle the meridian
        x_coast, y_coast, zcoast = ot.transform_points(ccrs.PlateCarree(), np.degrees(mlon), mlat).T
        plt.plot(x_coast, y_coast, zorder=0, color='grey', linewidth=0.5, alpha=0.6)

    return ax, ot, 'mag', fig
