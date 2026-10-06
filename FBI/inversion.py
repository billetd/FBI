"""
Cut down version of lompe's Emodel.run_inversion(), for line-of-sight convection data only.
Solves for just the model vector, skipping the posterior covariance and resolution matrices
which FBI doesn't use.
"""
import numpy as np
import scipy.linalg
from secsy import cubedsphere as cs


def prepare_inversion(model, l1=10, l2=0.1, perimeter_width=10):
    """
    Build the parts of the inversion that are the same for every scan
    :param model: lompe Emodel
    :param l1: float - Damping parameter for the model norm
    :param l2: float - Damping parameter for variation in the magnetic eastward direction
    :param perimeter_width: int - Cells added around grid_J for the data density grid.
                            Data outside of it are not used.
    :return: dict with the data density grid ('biggrid') and the regularisation matrix ('LTL')
    """

    # Expanded grid for calculation of data density, as in run_inversion()
    biggrid = cs.CSgrid(model.grid_J.projection,
                        model.grid_J.L + 2 * perimeter_width * model.grid_J.Lres,
                        model.grid_J.W + 2 * perimeter_width * model.grid_J.Wres,
                        model.grid_J.Lres, model.grid_J.Wres,
                        R=model.R)

    # Roughening matrix, as reg_E() in run_inversion() with l3 = 0
    LTL = 0
    if l1 > 0:
        LTL_l1 = np.eye(model.grid_E.size)
        LTL += l1 * LTL_l1 / np.median(LTL_l1.diagonal())
    if l2 > 0:
        LTL += l2 * model.LTLe / np.median(model.LTLe.diagonal())

    return {'biggrid': biggrid, 'LTL': LTL}


def solve_los(model, setup, lon, lat, vlos, le, ln, error):
    """
    Solve for the model vector given line-of-sight velocities
    :param model: lompe Emodel
    :param setup: dict from prepare_inversion()
    :param lon: array - Geographic longitudes of the data [degrees]
    :param lat: array - Geographic latitudes of the data [degrees]
    :param vlos: array - Line-of-sight speeds [m/s]
    :param le: array - Eastward components of the line-of-sight unit vectors
    :param ln: array - Northward components of the line-of-sight unit vectors
    :param error: array - Measurement errors [m/s]
    :return: (m, used) - The model vector, or None if there was too little data to fit,
             and a boolean array of which data points went into the fit
    """

    biggrid = setup['biggrid']

    # Same data selection as lompe: drop NaNs, then anything outside biggrid
    used = np.isfinite(vlos)
    used[used] = biggrid.ingrid(lon[used], lat[used])
    if used.sum() <= 1:
        return None, used

    lon, lat, vlos, le, ln, error = (x[used] for x in (lon, lat, vlos, le, ln, error))

    # Project the velocity matrices onto the line of sight
    Ge, Gn = model._v_matrix(lon=lon, lat=lat)
    G = Ge * le.reshape((-1, 1)) + Gn * ln.reshape((-1, 1))

    # Weights inversely proportional to data density
    bincount = biggrid.count(lon, lat)
    i, j = biggrid.bin_index(lon, lat)
    spatial_weight = 1. / np.maximum(bincount[i, j], 1)
    spatial_weight[i == -1] = 1
    w = spatial_weight * 1 / (error ** 2)

    GTG = G.T.dot(G * w.reshape((-1, 1)))
    GTd = G.T.dot(w * vlos)
    GG = GTG + setup['LTL'] * np.median(np.diagonal(GTG))

    try:
        c, lower = scipy.linalg.cho_factor(GG, lower=True)
        m = scipy.linalg.cho_solve((c, lower), GTd)
    except scipy.linalg.LinAlgError:
        m = scipy.linalg.lstsq(GG, GTd, cond=None, lapack_driver='gelsy')[0]

    return m, used
