"""
The values saved for each scan, and reading and writing them to FBI hdf5 files
"""
import numpy as np
import h5py
import datetime as dt

# How each dataset is saved, in the order they are written. Velocities and potentials are
# rounded to whole numbers (scaleoffset=0) to save space.
_ROUNDED = dict(compression="gzip", chunks=True, shuffle=True, scaleoffset=0, compression_opts=4)
_COMPRESSED = dict(compression="gzip")
_DATASETS = (
    # Fit vectors for the lompe grid
    ('v_e_model', _ROUNDED), ('v_n_model', _ROUNDED), ('mlats_model', _COMPRESSED), ('mlons_model', _COMPRESSED),
    # Line-of-sight data going into the fit
    ('v_e_los', _ROUNDED), ('v_n_los', _ROUNDED), ('mlats_los', _COMPRESSED), ('mlons_los', _COMPRESSED),
    ('rids', _COMPRESSED),
    # Fit vectors, at the locations of the SuperDARN equal area grid
    ('v_e_darngrid', _ROUNDED), ('v_n_darngrid', _ROUNDED), ('mlats_darngrid', _COMPRESSED),
    ('mlons_darngrid', _COMPRESSED),
    # Electric potential on the Lompe grid
    ('e_pot_model', _ROUNDED),
    # Coordinates of the boundary of the fit
    ('bound_mlats', _COMPRESSED), ('bound_mlons', _COMPRESSED),
)
_TIME_DATASETS = ('scan_year', 'scan_month', 'scan_day', 'scan_hour', 'scan_minute', 'scan_second',
                  'scan_millisec')


def output_geometry(model, apex, darn_grid_stuff):
    """
    Everything in lompe_extract() that depends only on the Lompe grid and the SuperDARN grid,
    not on the data of a particular scan
    :param model: FBI.inversion.Model
    :param apex: apexpy.Apex object
    :param darn_grid_stuff: dict from FBI.grid.sdarn_grid()
    :return: dict of arrays
    """

    # Restrict the superdarn grid points to those only in the lompe grid
    ingrid = model.grid_E.ingrid(darn_grid_stuff['glons_darngrid'], darn_grid_stuff['glats_darngrid'])
    glats_darngrid = darn_grid_stuff['glats_darngrid'][ingrid]
    glons_darngrid = darn_grid_stuff['glons_darngrid'][ingrid]

    # Model grid points
    glats_model, glons_model = model.grid_J.lat.flatten(), model.grid_J.lon.flatten()
    mlats_model, mlons_model = apex.geo2apex(glats_model, glons_model, 300)

    # Locations of the model boundary
    lon, lat = model.grid_J.lon, model.grid_J.lat
    bound_lons = np.concatenate((lon[0, :], lon[:, -1], np.flip(lon[-1, :]), np.flip(lon[:, 0])))
    bound_lats = np.concatenate((lat[0, :], lat[:, -1], np.flip(lat[-1, :]), np.flip(lat[:, 0])))
    bound_mlats, bound_mlons = apex.geo2apex(bound_lats, bound_lons, 300)

    return {'mlats_darngrid': darn_grid_stuff['mlats_darngrid'][ingrid],
            'mlons_darngrid': darn_grid_stuff['mlons_darngrid'][ingrid],
            'mlats_model': mlats_model, 'mlons_model': mlons_model,
            'bound_mlats': bound_mlats, 'bound_mlons': bound_mlons,
            # Apex base vectors
            'f_model': apex.basevectors_qd(glats_model, glons_model, 300, coords='geo'),
            'f_darngrid': apex.basevectors_qd(glats_darngrid, glons_darngrid, 300, coords='geo'),
            # SECS matrices. The velocities and potential are just these dotted with the model vector.
            'v_matrix_model': model.v_matrix(),
            'v_matrix_darngrid': model.v_matrix(glons_darngrid, glats_darngrid),
            'pot_matrix_model': model.potential_matrix()}


def _to_qd(f, v_e_geo, v_n_geo):
    """
    Rotate geographic velocity components into the magnetic (quasi-dipole) frame using
    apex base vectors. The direction is the motion across the QD grid, from the gradients of
    QD longitude and latitude (f2 x k and k x f1, F times g1 and g2 in Richmond (1995)
    equations (6.3) and (6.4)). The QD grid lines aren't perpendicular, so projecting onto
    f1 and f2 directly would skew the vectors. The speed is kept as the geographic speed.

    :param f: (f1, f2) tuple as returned by apexpy.Apex.basevectors_qd()
    :param v_e_geo: eastward velocity components
    :param v_n_geo: northward velocity components
    :return: (v_e_mag, v_n_mag)
    """
    f1, f2 = f
    v_e_mag = f2[1] * v_e_geo - f2[0] * v_n_geo
    v_n_mag = -f1[1] * v_e_geo + f1[0] * v_n_geo

    # Back to the geographic speed
    speed_geo, speed_mag = np.hypot(v_e_geo, v_n_geo), np.hypot(v_e_mag, v_n_mag)
    scale = np.divide(speed_geo, speed_mag, out=np.zeros_like(speed_mag), where=speed_mag != 0)
    return v_e_mag * scale, v_n_mag * scale


def lompe_extract(geometry, m, los, scan_time):
    """
    Potentials and velocities at good points for later plotting. These are the values saved
    to the hdf5 file.
    :param geometry: dict from output_geometry()
    :param m: model vector from FBI.inversion.Model.solve()
    :param los: dict of the data that went into the fit: geographic velocity 'v_e_geo' and 'v_n_geo',
                position 'mlats' and 'mlons', QD base vectors 'f' and station ids 'rids'
    :param scan_time: datetime of the scan
    :return: dict, one entry per dataset
    """

    # Model and darngrid velocities, rotated into the magnetic frame along with the data
    v_e_model, v_n_model = _to_qd(geometry['f_model'], *(V.dot(m) for V in geometry['v_matrix_model']))
    v_e_darngrid, v_n_darngrid = _to_qd(geometry['f_darngrid'], *(V.dot(m) for V in geometry['v_matrix_darngrid']))
    v_e_los, v_n_los = _to_qd(los['f'], los['v_e_geo'], los['v_n_geo'])

    return {'v_e_model': v_e_model, 'v_n_model': v_n_model,
            'mlats_model': geometry['mlats_model'], 'mlons_model': geometry['mlons_model'],
            'v_e_los': v_e_los, 'v_n_los': v_n_los, 'rids': los['rids'],
            'mlats_los': los['mlats'], 'mlons_los': los['mlons'],
            'v_e_darngrid': v_e_darngrid, 'v_n_darngrid': v_n_darngrid,
            'mlats_darngrid': geometry['mlats_darngrid'], 'mlons_darngrid': geometry['mlons_darngrid'],
            'e_pot_model': geometry['pot_matrix_model'].dot(m),
            'bound_mlats': geometry['bound_mlats'], 'bound_mlons': geometry['bound_mlons'],
            'scan_year': scan_time.year, 'scan_month': scan_time.month, 'scan_day': scan_time.day,
            'scan_hour': scan_time.hour, 'scan_minute': scan_time.minute,
            'scan_second': scan_time.second, 'scan_millisec': scan_time.microsecond}


def fbi_hdf5_name(timerange):
    """
    Name of the output file for a given timerange
    :param timerange: list[datetime] - Start and end times
    :return: str
    """

    return ('FBI_' + timerange[0].strftime("%Y%m%d%H%M%S") + '_'
            + timerange[1].strftime("%Y%m%d%H%M%S") + ".hdf5")


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

        grp = self._f.create_group(str(index))
        for name, options in _DATASETS:
            grp.create_dataset(name, data=lompe[name], **options)
        for name in _TIME_DATASETS:
            grp.create_dataset(name, shape=1, data=lompe[name])
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
    Load the data saved by FBIWriter or fbi_save_hdf5()
    Use timerange as a datetime tuple to only read in between two times
    :param file: Path to the FBI hdf5 file
    :param timerange: Optional. [start_time, end_time] datetime objects from a period of time to read in.
    :param as_arrays: Keep the datasets as numpy arrays rather than converting to lists
    :return: lompes: list of dictionaries containing the data
    """

    print('Reading: ' + file)

    lompes = []
    with h5py.File(file, "r") as f:

        # Groups are named by scan index, which hdf5 sorts as strings
        for name in sorted(f.keys(), key=int):
            group = f[name]

            # Check if it's within the timerange, if using
            if timerange:
                this_time = dt.datetime(*(group[key][0] for key in _TIME_DATASETS[:6]))
                if not timerange[0] <= this_time < timerange[1]:
                    continue

            this_record = {}
            for dataset in group:
                values = group[dataset][()]
                this_record[dataset] = values if as_arrays else values.tolist()
            lompes.append(this_record)

    return lompes
