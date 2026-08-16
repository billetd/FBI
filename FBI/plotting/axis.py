"""
Functions for building the map axes that FBI.plotting.plot draws onto.
"""
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
import numpy as np
import matplotlib.path as mpltPath
from matplotlib import pyplot as plt
from polplot import Polarplot
from shapely.geometry import MultiLineString


def get_local_axis(apex, center_lon=-61.15, center_lat=83, lon_bounds=(-82, 64), lat_bounds=(56, 71),
                   lat_gridline_range=(50, 91, 10), lon_gridline_range=(-180, 181, 45),
                   coastline_resolution='50m', refh=300):
    """
    Build an orthographic map axis in magnetic coordinates, centred over the Borealis Canada radars.
    :param apex: apexpy.Apex instance, used to convert coastlines from geographic to magnetic coordinates.
    :param center_lon: Longitude the orthographic projection is centred on.
    :param center_lat: Latitude the orthographic projection is centred on.
    :param lon_bounds: (min, max) longitude of the visible map window.
    :param lat_bounds: (min, max) latitude of the visible map window.
    :param lat_gridline_range: (start, stop, step) passed to np.arange() for latitude gridlines.
    :param lon_gridline_range: (start, stop, step) passed to np.arange() for longitude gridlines.
    :param coastline_resolution: Natural Earth coastline resolution to draw.
    :param refh: Apex reference height in km, used for the coastline coordinate conversion.
    :return: ax, projection, coord ('mag'), fig
    """

    # Set the projection to orthographic
    projection = ccrs.Orthographic(center_lon, center_lat)
    fig, ax = plt.subplots(subplot_kw={'projection': projection}, gridspec_kw={'wspace': 0.05, 'hspace': 0.05})

    # Set up the plot
    xs, ys, _ = projection.transform_points(
        ccrs.PlateCarree(),
        np.array((lon_bounds[0], lon_bounds[1])),
        np.array((lat_bounds[0], lat_bounds[1])),
    ).T

    # For checking if coastline points are in the plot window later
    window = mpltPath.Path([[xs[0], ys[0]], [xs[0], ys[1]], [xs[1], ys[1]], [xs[1], ys[0]]])

    # Setup the axis to look nice
    ax.set_xlim(xs)
    ax.set_ylim(ys)
    gl = ax.gridlines(draw_labels=True, y_inline=True, crs=ccrs.PlateCarree(), zorder=5, color='grey')
    gl.ylocator = mticker.FixedLocator(np.arange(*lat_gridline_range))
    gl.ylabel_style = {'color': 'black'}
    gl.xlocator = mticker.FixedLocator(np.arange(*lon_gridline_range))
    gl.top_labels = False
    gl.bottom_labels = False
    gl.xlines = False

    # Read in the geometry object of the coastlines
    coastlines = cfeature.NaturalEarthFeature('physical', 'coastline', coastline_resolution,
                                              color='k', zorder=2.0)

    # Plot coastlines
    for shape in coastlines.geometries():
        if isinstance(shape, MultiLineString):  # Don't plot multi geoms as it breaks
            continue
        glons, glats = shape.coords.xy[0], shape.coords.xy[1]
        mlats, mlons = apex.geo2apex(glats, glons, refh)
        x_coast, y_coast, _ = projection.transform_points(ccrs.PlateCarree(), mlons, mlats).T
        if any(window.contains_points(list(zip(x_coast, y_coast)))):  # Only plot if visible in the window
            ax.plot(x_coast, y_coast, zorder=0, color='grey', linewidth=0.5, alpha=0.6)

    return ax, projection, 'mag', fig


def get_polar_axis(time, apex, minlat=50, linewidth=0.7, coastline_linewidth=0.5, coastline_resolution='50m'):
    """
    Build a full magnetic-local-time polar map axis.
    :param time: datetime of the scan, used to position coastlines correctly.
    :param apex: apexpy.Apex instance, used by Polarplot to draw coastlines.
    :param minlat: Lowest latitude shown on the polar plot.
    :param linewidth: Line width for the plot's gridlines and low-latitude boundary.
    :param coastline_linewidth: Line width for coastlines.
    :param coastline_resolution: Natural Earth coastline resolution to draw.
    :return: pax (Polarplot), coord ('mlt'), fig
    """

    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

    pax = Polarplot(ax, minlat=minlat, linewidth=linewidth)
    pax.coastlines(time=time, mag=apex, linewidth=coastline_linewidth, resolution=coastline_resolution,
                  color='grey', alpha=0.5)

    lowlat_mlts = np.linspace(0, 24, num=360)
    lowlat_lats = np.full(360, minlat)

    # Polarplot's default low-latitude boundary is dashed; draw a solid one instead
    pax.plot(lowlat_lats, lowlat_mlts, color='black')

    return pax, 'mlt', fig
