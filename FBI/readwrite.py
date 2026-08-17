import numpy as np
import h5py
import datetime as dt
from secsy import get_SECS_J_G_matrices

# Cache of the geometry that does not change between scans, built once and then inherited
# by the workers. See _build_geometry_cache() for what goes in it.
_geometry_cache = None


def prime_geometry_cache(model, apex, darn_grid_stuff):
    """
    Build the static geometry up front, so the workers inherit it rather than each
    building its own. Everything it needs is on the Emodel before any inversion is run.
    :param model: lompe Emodel, used only for its grids
    :param apex: apexpy.Apex object
    :param darn_grid_stuff: dict from FBI.grid.sdarn_grid()
    """

    global _geometry_cache
    _geometry_cache = _build_geometry_cache(model, apex, darn_grid_stuff)


def reset_geometry_cache():
    """
    Drop the cached geometry. It isn't keyed on anything, so it has to go whenever the
    grid or the apex epoch changes.
    """

    global _geometry_cache
    _geometry_cache = None


def _build_geometry_cache(scan_lompe, apex, darn_grid_stuff):
    """
    Precompute everything in lompe_extract() that depends only on the Lompe grid and the
    SuperDARN grid, not on the data of a particular scan.

    :param scan_lompe: lompe Emodel, used only for its grids
    :param apex: apexpy.Apex object
    :param darn_grid_stuff: dict from FBI.grid.sdarn_grid()
    :return: dict of cached arrays
    """

    # Restrict the superdarn grid points to those only in the lompe grid
    ingrid = scan_lompe.grid_E.ingrid(darn_grid_stuff['glons_darngrid'], darn_grid_stuff['glats_darngrid'])
    glats_darngrid = darn_grid_stuff['glats_darngrid'][ingrid]
    glons_darngrid = darn_grid_stuff['glons_darngrid'][ingrid]

    # Model grid points
    glats_model, glons_model = scan_lompe.grid_J.lat.flatten(), scan_lompe.grid_J.lon.flatten()

    # Locations of the model boundary
    bound_lons = np.concatenate((scan_lompe.grid_J.lon[0, :], scan_lompe.grid_J.lon[:, -1],
                                 np.flip(scan_lompe.grid_J.lon[-1, :]), np.flip(scan_lompe.grid_J.lon[:, 0])))
    bound_lats = np.concatenate((scan_lompe.grid_J.lat[0, :], scan_lompe.grid_J.lat[:, -1],
                                 np.flip(scan_lompe.grid_J.lat[-1, :]), np.flip(scan_lompe.grid_J.lat[:, 0])))
    bound_mlats, bound_mlons = apex.geo2apex(bound_lats, bound_lons, 300)

    mlats_model, mlons_model = apex.geo2apex(glats_model, glons_model, 300)

    cache = {
        'ingrid': ingrid,
        'glats_darngrid': glats_darngrid,
        'glons_darngrid': glons_darngrid,
        'mlats_darngrid': darn_grid_stuff['mlats_darngrid'][ingrid],
        'mlons_darngrid': darn_grid_stuff['mlons_darngrid'][ingrid],
        'mlats_model': mlats_model,
        'mlons_model': mlons_model,
        'bound_mlats': bound_mlats,
        'bound_mlons': bound_mlons,
        # Apex base vectors for the two static point sets
        'f_model': apex.basevectors_qd(glats_model, glons_model, 300, coords='geo'),
        'f_darngrid': apex.basevectors_qd(glats_darngrid, glons_darngrid, 300, coords='geo'),
        # SECS matrices. v() and E_pot() are just these dotted with the model vector.
        'v_matrix_model': scan_lompe._v_matrix(),
        'v_matrix_darngrid': scan_lompe._v_matrix(lon=glons_darngrid, lat=glats_darngrid),
        'pot_matrix_model': get_SECS_J_G_matrices(scan_lompe.lat_J, scan_lompe.lon_J,
                                                  scan_lompe.lat_E, scan_lompe.lon_E,
                                                  current_type='potential',
                                                  RI=scan_lompe.R,
                                                  singularity_limit=scan_lompe.secs_singularity_limit),
    }

    return cache


def _to_qd(f, v_e_geo, v_n_geo):
    """
    Rotate geographic velocity components into the magnetic (quasi-dipole) frame using
    apex base vectors. Richmond (1995) equations (7.12) and (7.13), but for velocities.

    :param f: (f1, f2) tuple as returned by apexpy.Apex.basevectors_qd()
    :param v_e_geo: eastward velocity components
    :param v_n_geo: northward velocity components
    :return: (v_e_mag, v_n_mag)
    """
    f1, f2 = f
    return (f1[0] * v_e_geo + f1[1] * v_n_geo,
            f2[0] * v_e_geo + f2[1] * v_n_geo)


def lompe_extract(scan_lompe, apex, scan_time, darn_grid_stuff, rids, use_cache=True):
    """
    Code to extract potentials and velocities at good points for later plotting
    These are the values saved to HDF5 later in fbi_write_hdf5()
    :param scan_lompe:
    :param apex:
    :param scan_time:
    :param darn_grid_stuff:
    :param rids: NOTE - I think this needs to be fixed to remove radars outside of lompe grid area
    :param use_cache: Reuse the static geometry between scans. Much faster, but holds a few
                      hundred MB of SECS matrices per process. Set False if RAM-limited.
    :return:
    """

    global _geometry_cache

    if use_cache:
        if _geometry_cache is None:
            _geometry_cache = _build_geometry_cache(scan_lompe, apex, darn_grid_stuff)
        geom = _geometry_cache
    else:
        geom = _build_geometry_cache(scan_lompe, apex, darn_grid_stuff)

    m = scan_lompe.m

    # Darngrid velocities
    Ve, Vn = geom['v_matrix_darngrid']
    v_e_geo_darngrid, v_n_geo_darngrid = Ve.dot(m), Vn.dot(m)

    # Model velocities
    Ve, Vn = geom['v_matrix_model']
    v_e_geo_model, v_n_geo_model = Ve.dot(m), Vn.dot(m)

    # Electric potential
    e_pot_model = geom['pot_matrix_model'].dot(m)

    # Data velocities and points. These are the only points that move between scans.
    v_e_geo_los, v_n_geo_los = (scan_lompe.data['convection'][0].values * scan_lompe.data['convection'][0].los_mag[0],
                                scan_lompe.data['convection'][0].values * scan_lompe.data['convection'][0].los_mag[1])
    glons_los, glats_los = (scan_lompe.data['convection'][0].coords['lon'],
                            scan_lompe.data['convection'][0].coords['lat'])
    mlats_los, mlons_los = apex.geo2apex(glats_los, glons_los, 300)

    # Rotate all three sets of velocities into the magnetic frame
    v_e_model, v_n_model = _to_qd(geom['f_model'], v_e_geo_model, v_n_geo_model)
    v_e_darngrid, v_n_darngrid = _to_qd(geom['f_darngrid'], v_e_geo_darngrid, v_n_geo_darngrid)
    v_e_los, v_n_los = _to_qd(apex.basevectors_qd(glats_los, glons_los, 300, coords='geo'),
                              v_e_geo_los, v_n_geo_los)

    data = {'v_e_model': v_e_model, 'v_n_model': v_n_model,
            'mlats_model': geom['mlats_model'], 'mlons_model': geom['mlons_model'],
            'v_e_los': v_e_los, 'v_n_los': v_n_los, 'rids': rids,
            'mlats_los': mlats_los, 'mlons_los': mlons_los,
            'v_e_darngrid': v_e_darngrid, 'v_n_darngrid': v_n_darngrid,
            'mlats_darngrid': geom['mlats_darngrid'], 'mlons_darngrid': geom['mlons_darngrid'],
            'e_pot_model': e_pot_model,
            'bound_mlats': geom['bound_mlats'], 'bound_mlons': geom['bound_mlons'],
            'scan_year': scan_time.year, 'scan_month': scan_time.month, 'scan_day': scan_time.day,
            'scan_hour': scan_time.hour, 'scan_minute': scan_time.minute,
            'scan_second': scan_time.second, 'scan_millisec': scan_time.microsecond}

    return data


def fbi_hdf5_name(timerange):
    """
    Name of the output file for a given timerange
    :param timerange: list[datetime] - Start and end times
    :return: str
    """

    return ('FBI_' + timerange[0].strftime("%Y%m%d%H%M%S") + '_'
            + timerange[1].strftime("%Y%m%d%H%M%S") + ".hdf5")


def _write_scan_group(f, counter, lompe):
    """
    Write one scan's output as a group of an open hdf5 file
    :param f: open h5py.File
    :param counter: int - scan index, used as the group name
    :param lompe: dict from lompe_extract()
    """

    grp = f.create_group(str(counter))

    # Fit vectors for the lompe grid
    grp.create_dataset("v_e_model", shape=(len(lompe['v_e_model'])), data=lompe['v_e_model'],
                       compression="gzip", chunks=True, shuffle=True, scaleoffset=0, compression_opts=4)
    grp.create_dataset("v_n_model", shape=(len(lompe['v_n_model'])), data=lompe['v_n_model'],
                       compression="gzip", chunks=True, shuffle=True, scaleoffset=0, compression_opts=4)
    grp.create_dataset("mlats_model", shape=(len(lompe['mlats_model'])), data=lompe['mlats_model'],
                       compression="gzip")
    grp.create_dataset("mlons_model", shape=(len(lompe['mlons_model'])), data=lompe['mlons_model'],
                       compression="gzip")

    # Line-of-sight data going into the fit
    grp.create_dataset("v_e_los", shape=(len(lompe['v_e_los'])), data=lompe['v_e_los'],
                       compression="gzip", chunks=True, shuffle=True, scaleoffset=0, compression_opts=4)
    grp.create_dataset("v_n_los", shape=(len(lompe['v_n_los'])), data=lompe['v_n_los'],
                       compression="gzip", chunks=True, shuffle=True, scaleoffset=0, compression_opts=4)
    grp.create_dataset("mlats_los", shape=(len(lompe['mlats_los'])), data=lompe['mlats_los'],
                       compression="gzip")
    grp.create_dataset("mlons_los", shape=(len(lompe['mlons_los'])), data=lompe['mlons_los'],
                       compression="gzip")
    grp.create_dataset("rids", shape=(len(lompe['rids'])), data=lompe['rids'], compression="gzip")

    # Fit vectors, at the locations of the SuperDARN equal area grid
    grp.create_dataset("v_e_darngrid", shape=(len(lompe['v_e_darngrid'])), data=lompe['v_e_darngrid'],
                       compression="gzip", chunks=True, shuffle=True, scaleoffset=0, compression_opts=4)
    grp.create_dataset("v_n_darngrid", shape=(len(lompe['v_n_darngrid'])), data=lompe['v_n_darngrid'],
                       compression="gzip", chunks=True, shuffle=True, scaleoffset=0, compression_opts=4)
    grp.create_dataset("mlats_darngrid", shape=(len(lompe['mlats_darngrid'])), data=lompe['mlats_darngrid'],
                       compression="gzip")
    grp.create_dataset("mlons_darngrid", shape=(len(lompe['mlons_darngrid'])), data=lompe['mlons_darngrid'],
                       compression="gzip")

    # Electric potential on the Lompe grid
    grp.create_dataset("e_pot_model", shape=(len(lompe['e_pot_model'])), data=lompe['e_pot_model'],
                       compression="gzip", chunks=True, shuffle=True, scaleoffset=0, compression_opts=4)

    # Coordinates of the boundary of the fit
    grp.create_dataset("bound_mlats", shape=(len(lompe['bound_mlats'])), data=lompe['bound_mlats'],
                       compression="gzip")
    grp.create_dataset("bound_mlons", shape=(len(lompe['bound_mlons'])), data=lompe['bound_mlons'],
                       compression="gzip")

    # Time info
    grp.create_dataset("scan_year", shape=1, data=lompe['scan_year'])
    grp.create_dataset("scan_month", shape=1, data=lompe['scan_month'])
    grp.create_dataset("scan_day", shape=1, data=lompe['scan_day'])
    grp.create_dataset("scan_hour", shape=1, data=lompe['scan_hour'])
    grp.create_dataset("scan_minute", shape=1, data=lompe['scan_minute'])
    grp.create_dataset("scan_second", shape=1, data=lompe['scan_second'])
    grp.create_dataset("scan_millisec", shape=1, data=lompe['scan_millisec'])


class FBIWriter:
    """
    Writes each scan to the hdf5 file as it is fitted, so a whole run of them never has
    to be held in memory at once. One group per scan index, skipping scans with no fit.

    Usage:
        with FBIWriter(timerange, lompe_dir) as writer:
            writer.write(index, lompe_data)
    """

    def __init__(self, timerange, lompe_dir):
        """
        :param timerange: list[datetime] - Start and end times, used for the file name
        :param lompe_dir: str - Directory to save the FBI output file in
        """

        if not lompe_dir.endswith('/'):
            lompe_dir += '/'
        self.path = lompe_dir + fbi_hdf5_name(timerange)
        self._f = None
        self.n_written = 0

    def __enter__(self):
        print('Writing to ' + self.path)
        self._f = h5py.File(self.path, "w")
        return self

    def write(self, index, lompe):
        """
        :param index: int - scan index, used as the group name
        :param lompe: dict from lompe_extract(), or None if the scan produced no fit
        """

        if lompe is None:
            return
        _write_scan_group(self._f, index, lompe)
        self.n_written += 1

    def __exit__(self, *exc):
        self._f.close()
        self._f = None
        return False


def fbi_save_hdf5(lompes, timerange, lompe_dir):
    """
    Save a complete list of scans into a hdf5 file, all at once
    :param lompes: list of lompe_extract() dicts, None where no fit was made
    :param timerange: list[datetime] - Start and end times
    :param lompe_dir: str - Directory to save the FBI output file in
    """

    with FBIWriter(timerange, lompe_dir) as writer:
        for counter, lompe in enumerate(lompes):
            writer.write(counter, lompe)


def fbi_load_hdf5(file, timerange=None, as_arrays=False):
    """
    Load the data saved by fbi_save_hdf5()
    Use timerange as a datetime tuple to only read in between two times
    :param file: Path to the FBI hdf5 file
    :param timerange: Optional. [start_time, end_time] datetime objects from a period of time to read in.
    :param as_arrays: Keep the datasets as numpy arrays rather than converting to lists
    :return: lompes: list of dictionaries containing the data
    """

    print('Reading: ' + file)

    with h5py.File(file, "r") as f:
        lompes = []
        groups_keys = list(f.keys())

        # Fix the wonky hdf5 sorting
        group_ints = [int(group) for group in groups_keys]
        group_ints.sort()
        groups = [str(group) for group in group_ints]

        # Iterate over records
        for group in groups:

            datasets = f[group].keys()

            # Check if it's within the timerange, if using
            if timerange:
                this_time = dt.datetime(f[group + '/' + 'scan_year'][0],
                                        f[group + '/' + 'scan_month'][0],
                                        f[group + '/' + 'scan_day'][0],
                                        f[group + '/' + 'scan_hour'][0],
                                        f[group + '/' + 'scan_minute'][0],
                                        f[group + '/' + 'scan_second'][0])
                if timerange[0] <= this_time < timerange[1]:
                    pass
                else:
                    continue

            this_record = {}
            # Iterate over keys
            for dataset in datasets:

                values = f[group + '/' + dataset][()]
                this_record[dataset] = values if as_arrays else values.tolist()

            # Append to list of dictionaries
            lompes.append(this_record)

    return lompes



