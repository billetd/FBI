import numpy as np
import h5py


def lompe_extract(scan_lompe, apex, scan_time, darn_grid_stuff):
    """
    Code to extract potentials and velocities at good points for later plotting
    These are the values saved to HDF5 later in fbi_write_hdf5()
    :param scan_lompe:
    :param apex:
    :param scan_time:
    :param darn_grid_stuff:
    :return:
    """

    # Restrict the superdarn grid points to those only in the lompe grid
    ingrid = scan_lompe.grid_E.ingrid(darn_grid_stuff['glons_darngrid'], darn_grid_stuff['glats_darngrid'])
    glats_darngrid = darn_grid_stuff['glats_darngrid'][ingrid]
    glons_darngrid = darn_grid_stuff['glons_darngrid'][ingrid]
    mlats_darngrid = darn_grid_stuff['mlats_darngrid'][ingrid]
    mlons_darngrid = darn_grid_stuff['mlons_darngrid'][ingrid]

    # Darngrid velocities
    v_e_geo_darngrid, v_n_geo_darngrid = scan_lompe.v(lon=glons_darngrid, lat=glats_darngrid)

    # Model velocities and points
    v_e_geo_model, v_n_geo_model = scan_lompe.v()
    glats_model, glons_model = scan_lompe.grid_J.lat.flatten(), scan_lompe.grid_J.lon.flatten()
    mlats_model, mlons_model = apex.geo2apex(glats_model, glons_model, 300)

    # Data velocities and points
    v_e_geo_los, v_n_geo_los = (scan_lompe.data['convection'][0].values * scan_lompe.data['convection'][0].los_mag[0],
                                scan_lompe.data['convection'][0].values * scan_lompe.data['convection'][0].los_mag[1])
    glons_los, glats_los = (scan_lompe.data['convection'][0].coords['lon'],
                            scan_lompe.data['convection'][0].coords['lat'])
    mlats_los, mlons_los = apex.geo2apex(glats_los, glons_los, 300)

    # Locations of the model boundary
    bound_lons = np.concatenate((scan_lompe.grid_J.lon[0, :], scan_lompe.grid_J.lon[:, -1],
                                 np.flip(scan_lompe.grid_J.lon[-1, :]), np.flip(scan_lompe.grid_J.lon[:, 0])))
    bound_lats = np.concatenate((scan_lompe.grid_J.lat[0, :], scan_lompe.grid_J.lat[:, -1],
                                 np.flip(scan_lompe.grid_J.lat[-1, :]), np.flip(scan_lompe.grid_J.lat[:, 0])))
    bound_mlats, bound_mlons = apex.geo2apex(bound_lats, bound_lons, 300)

    # Electric potential
    e_pot_model = scan_lompe.E_pot()

    # Time
    hour = scan_time.hour
    minute = scan_time.minute
    second = scan_time.second
    year = scan_time.year
    month = scan_time.month
    day = scan_time.day
    millisec = scan_time.microsecond

    # Convert model velocities to mag frame
    # Get apex base vectors in geographic
    f1, f2 = apex.basevectors_qd(glats_model, glons_model, 300, coords='geo')
    # Rotate the geo veolocity vectors into magnetic using the base vectors
    # Richmond (1995) equations (7.12) and (7.13) but for velocity vectors
    v_geo_grid = np.vstack((v_e_geo_model, v_n_geo_model))
    v_e_model = np.einsum('ij,ij->j', f1, v_geo_grid)
    v_n_model = np.einsum('ij,ij->j', f2, v_geo_grid)

    # Convert data velocities to mag frame
    # Get apex base vectors in geographic
    f1, f2 = apex.basevectors_qd(glats_los, glons_los, 300, coords='geo')
    # Rotate the geo veolocity vectors into magnetic using the base vectors
    # Richmond (1995) equations (7.12) and (7.13) but for velocity vectors
    v_geo_grid = np.vstack((v_e_geo_los, v_n_geo_los))
    v_e_los = np.einsum('ij,ij->j', f1, v_geo_grid)
    v_n_los = np.einsum('ij,ij->j', f2, v_geo_grid)

    # darngrid velocities to mag frame
    f1, f2 = apex.basevectors_qd(glats_darngrid, glons_darngrid, 300, coords='geo')
    # Rotate the geo veolocity vectors into magnetic using the base vectors
    # Richmond (1995) equations (7.12) and (7.13) but for velocity vectors
    v_geo_grid = np.vstack((v_e_geo_darngrid, v_n_geo_darngrid))
    v_e_darngrid = np.einsum('ij,ij->j', f1, v_geo_grid)
    v_n_darngrid = np.einsum('ij,ij->j', f2, v_geo_grid)

    data = {'v_e_model': v_e_model.tolist(), 'v_n_model': v_n_model.tolist(),
            'mlats_model': mlats_model.tolist(), 'mlons_model': mlons_model.tolist(),
            'v_e_los': v_e_los.tolist(), 'v_n_los': v_n_los.tolist(),
            'mlats_los': mlats_los.tolist(), 'mlons_los': mlons_los.tolist(),
            'v_e_darngrid': v_e_darngrid.tolist(), 'v_n_darngrid': v_n_darngrid.tolist(),
            'mlats_darngrid': mlats_darngrid.tolist(), 'mlons_darngrid': mlons_darngrid.tolist(),
            'e_pot_model': e_pot_model.tolist(),
            'bound_mlats': bound_mlats.tolist(), 'bound_mlons': bound_mlons.tolist(),
            'scan_year': year, 'scan_month': month, 'scan_day': day, 'scan_hour': hour,
            'scan_minute': minute, 'scan_second': second, 'scan_millisec': millisec}

    return data


def fbi_save_hdf5(lompes, timerange, lompe_dir):
    """
    Save the output from process() into a hdf5 file
    :param lompes:
    :param timerange:
    :param lompe_dir:
    :return:
    """

    # TODO: Records don't appear to be in time order. Think is is because of using ray, or actually,
    # maybe because it's doing it in integer order (e.g. 12, 22, etc)
    # If plotting everything it all comes out eventually, but probably should be fixed

    # Dump data to hdf5 file
    print('Writing to file...')
    hdf5name = 'FBI_' + timerange[0].strftime("%Y%m%d%H%M%S") + '_' + timerange[1].strftime("%Y%m%d%H%M%S") + ".hdf5"
    with h5py.File(lompe_dir + hdf5name, "w") as f:
        for counter, lompe in enumerate(lompes):
            if lompe is not None:
                grp = f.create_group(str(counter))
                grp.create_dataset("v_e_model", shape=(len(lompe['v_e_model'])), data=lompe['v_e_model'],
                                   compression="gzip", chunks=True, shuffle=True, scaleoffset=0, compression_opts=9)
                grp.create_dataset("v_n_model", shape=(len(lompe['v_n_model'])), data=lompe['v_n_model'],
                                   compression="gzip", chunks=True, shuffle=True, scaleoffset=0, compression_opts=9)
                grp.create_dataset("mlats_model", shape=(len(lompe['mlats_model'])), data=lompe['mlats_model'],
                                   compression="gzip", chunks=True, shuffle=True, scaleoffset=0, compression_opts=9)
                grp.create_dataset("mlons_model", shape=(len(lompe['mlons_model'])), data=lompe['mlons_model'],
                                   compression="gzip", chunks=True, shuffle=True, scaleoffset=0, compression_opts=9)

                grp.create_dataset("v_e_los", shape=(len(lompe['v_e_los'])), data=lompe['v_e_los'],
                                   compression="gzip", chunks=True, shuffle=True, scaleoffset=0, compression_opts=9)
                grp.create_dataset("v_n_los", shape=(len(lompe['v_n_los'])), data=lompe['v_n_los'],
                                   compression="gzip", chunks=True, shuffle=True, scaleoffset=0, compression_opts=9)
                grp.create_dataset("mlats_los", shape=(len(lompe['mlats_los'])), data=lompe['mlats_los'],
                                   compression="gzip", chunks=True, shuffle=True, scaleoffset=0, compression_opts=9)
                grp.create_dataset("mlons_los", shape=(len(lompe['mlons_los'])), data=lompe['mlons_los'],
                                   compression="gzip", chunks=True, shuffle=True, scaleoffset=0, compression_opts=9)

                grp.create_dataset("v_e_darngrid", shape=(len(lompe['v_e_darngrid'])), data=lompe['v_e_darngrid'],
                                   compression="gzip", chunks=True, shuffle=True, scaleoffset=0, compression_opts=9)
                grp.create_dataset("v_n_darngrid", shape=(len(lompe['v_n_darngrid'])), data=lompe['v_n_darngrid'],
                                   compression="gzip", chunks=True, shuffle=True, scaleoffset=0, compression_opts=9)
                grp.create_dataset("mlats_darngrid", shape=(len(lompe['mlats_darngrid'])), data=lompe['mlats_darngrid'],
                                   compression="gzip", chunks=True, shuffle=True, scaleoffset=0, compression_opts=9)
                grp.create_dataset("mlons_darngrid", shape=(len(lompe['mlons_darngrid'])), data=lompe['mlons_darngrid'],
                                   compression="gzip", chunks=True, shuffle=True, scaleoffset=0, compression_opts=9)

                grp.create_dataset("e_pot_model", shape=(len(lompe['e_pot_model'])), data=lompe['e_pot_model'],
                                   compression="gzip", chunks=True, shuffle=True, scaleoffset=0, compression_opts=9)

                grp.create_dataset("bound_mlats", shape=(len(lompe['bound_mlats'])), data=lompe['bound_mlats'],
                                   compression="gzip", chunks=True, shuffle=True, scaleoffset=0, compression_opts=9)

                grp.create_dataset("scan_year", shape=1, data=lompe['scan_year'], compression="gzip",
                                   chunks=True, shuffle=True, scaleoffset=0, compression_opts=9)
                grp.create_dataset("scan_month", shape=1, data=lompe['scan_month'], compression="gzip",
                                   chunks=True, shuffle=True, scaleoffset=0, compression_opts=9)
                grp.create_dataset("scan_day", shape=1, data=lompe['scan_day'], compression="gzip",
                                   chunks=True, shuffle=True, scaleoffset=0, compression_opts=9)
                grp.create_dataset("scan_hour", shape=1, data=lompe['scan_hour'], compression="gzip",
                                   chunks=True, shuffle=True, scaleoffset=0, compression_opts=9)
                grp.create_dataset("scan_minute", shape=1, data=lompe['scan_minute'],
                                   compression="gzip", chunks=True, shuffle=True, scaleoffset=0, compression_opts=9)
                grp.create_dataset("scan_second", shape=1, data=lompe['scan_second'],
                                   compression="gzip", chunks=True, shuffle=True, scaleoffset=0, compression_opts=9)
                grp.create_dataset("scan_millisec", shape=1, data=lompe['scan_millisec'],
                                   compression="gzip", chunks=True, shuffle=True, scaleoffset=0, compression_opts=9)


def fbi_load_hdf5(file):
    """
    Load the data saved by fbi_save_hdf5()
    :param file:
    :return:
    """

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
            this_record = {}

            # Iterate over keys
            for dataset in datasets:

                this_record[dataset] = f[group + '/' + dataset][()].tolist()

            # Append to list of dictionaries
            lompes.append(this_record)

    return lompes

