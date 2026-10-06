"""
Cut down versions of lompe's Emodel and Emodel.run_inversion(), for line-of-sight convection data only.
Builds just the matrices FBI uses, and solves for just the model vector, skipping the posterior
covariance and resolution matrices.
"""
import datetime as dt
import apexpy
import numpy as np
import scipy.linalg
from ppigrf import igrf
from secsy import cubedsphere as cs
from secsy import get_SECS_J_G_matrices

RE = 6371.2e3  # Earth radius, as in lompe [m]


class Model:

    def __init__(self, grid, l1=10, l2=0.1, ew_regularization_limit=(50, 75), perimeter_width=10,
                 epoch=dt.datetime(2015, 1, 1)):
        """
        :param grid: secsy CSgrid - The lompe grid (grid_J in lompe)
        :param l1: float - Damping parameter for the model norm
        :param l2: float - Damping parameter for variation in the magnetic eastward direction
        :param ew_regularization_limit: (float, float) - Magnetic latitudes between which the east-west
                                        regularisation is reduced to zero towards the pole
        :param perimeter_width: int - Cells added around grid_J for the data density grid.
                                Data outside of it are not used.
        :param epoch: datetime - For the main field and the magnetic east direction. lompe's default.
        """

        # Inner and outer grids, as in Emodel
        self.grid_J = grid
        self.R = grid.R
        xi_e = np.hstack((grid.xi_mesh[0], grid.xi_mesh[0, -1] + grid.dxi)) - grid.dxi / 2
        eta_e = np.hstack((grid.eta_mesh[:, 0], grid.eta_mesh[-1, 0] + grid.deta)) - grid.deta / 2
        self.grid_E = cs.CSgrid(cs.CSprojection(grid.projection.position, grid.projection.orientation),
                                grid.L + grid.Lres, grid.W + grid.Wres, grid.Lres, grid.Wres,
                                edges=(xi_e, eta_e), R=self.R)
        self.lat_J, self.lon_J = np.ravel(self.grid_J.lat), np.ravel(self.grid_J.lon)
        self.lat_E, self.lon_E = np.ravel(self.grid_E.lat), np.ravel(self.grid_E.lon)
        self.secs_singularity_limit = np.min([grid.Wres, grid.Lres]) / 2

        # Main field on grid_E [T]
        refh = (self.R - RE) * 1e-3
        Be, Bn, Bu = igrf(self.lon_E, self.lat_E, refh, epoch)
        Be, Bn, Bu = Be * 1e-9, Bn * 1e-9, Bu * 1e-9
        self.B0 = np.sqrt(Be ** 2 + Bn ** 2 + Bu ** 2).reshape((1, -1))
        self.Bu = Bu.reshape((1, -1))
        hemisphere = -np.sign(self.Bu.flatten()[0])

        # Derivative in the magnetic eastward direction, as in Emodel.compute_L_matrices().
        # This leaves apexpy at the lompe epoch, so make the Apex for the data afterwards.
        De, Dn = self.grid_E.get_Le_Ln()
        apx = apexpy.Apex(epoch, refh=refh)
        mlat, _ = apx.geo2apex(self.lat_E, self.lon_E, refh)
        f1, _ = apx.basevectors_qd(self.lat_E, self.lon_E, refh)
        f1 = (f1 / np.linalg.norm(f1, axis=0)).reshape((-1, 1))
        a, b = ew_regularization_limit
        lat = hemisphere * mlat
        lat_w = np.where(lat < a, 1, np.where(lat > b, 0, (b - lat) / (b - a))).reshape((-1, 1))
        Le = (De * f1[0] + Dn * f1[1]) * lat_w
        LTLe = Le.T.dot(Le)

        # Regularisation, as reg_E() in run_inversion() with l3 = 0
        self.LTL = 0
        if l1 > 0:
            LTL_l1 = np.eye(self.grid_E.size)
            self.LTL += l1 * LTL_l1 / np.median(LTL_l1.diagonal())
        if l2 > 0:
            self.LTL += l2 * LTLe / np.median(LTLe.diagonal())

        # Expanded grid for calculation of data density, as in run_inversion()
        self.biggrid = cs.CSgrid(grid.projection,
                                 grid.L + 2 * perimeter_width * grid.Lres,
                                 grid.W + 2 * perimeter_width * grid.Wres,
                                 grid.Lres, grid.Wres, R=self.R)

    def v_matrix(self, lon=None, lat=None):
        """
        Matrices that give the eastward and northward velocities when dotted with the model vector,
        as Emodel._v_matrix()
        :param lon: array - Geographic longitudes [degrees]. Default is grid_J
        :param lat: array - Geographic latitudes [degrees]. Default is grid_J
        :return: (Ve, Vn)
        """

        if lon is None:
            lon, lat = self.lon_J, self.lat_J

        Ee, En = get_SECS_J_G_matrices(lat, lon, self.lat_E, self.lon_E, current_type='curl_free',
                                       RI=self.R, singularity_limit=self.secs_singularity_limit)
        return En * self.Bu / self.B0 ** 2, -Ee * self.Bu / self.B0 ** 2

    def potential_matrix(self):
        """
        Matrix that gives the electric potential on grid_J when dotted with the model vector
        """

        return get_SECS_J_G_matrices(self.lat_J, self.lon_J, self.lat_E, self.lon_E, current_type='potential',
                                     RI=self.R, singularity_limit=self.secs_singularity_limit)

    def los_matrix(self, lon, lat, le, ln):
        """
        Matrix that gives the velocity along the line of sight when dotted with the model vector
        :param lon: array - Geographic longitudes [degrees]
        :param lat: array - Geographic latitudes [degrees]
        :param le: array - Eastward components of the line-of-sight unit vectors
        :param ln: array - Northward components of the line-of-sight unit vectors
        :return: array, one row per point
        """

        Ve, Vn = self.v_matrix(lon, lat)
        return Ve * le.reshape((-1, 1)) + Vn * ln.reshape((-1, 1))

    def solve(self, G, lon, lat, d, error):
        """
        Solve for the model vector, with data weighted inversely by their density
        :param G: array - Data kernel, one row per point. See los_matrix()
        :param lon: array - Geographic longitudes of the data [degrees]
        :param lat: array - Geographic latitudes of the data [degrees]
        :param d: array - The data [m/s]
        :param error: array - Measurement errors [m/s]
        :return: The model vector
        """

        bincount = self.biggrid.count(lon, lat)
        i, j = self.biggrid.bin_index(lon, lat)
        spatial_weight = 1. / np.maximum(bincount[i, j], 1)
        spatial_weight[i == -1] = 1
        w = spatial_weight * 1 / (error ** 2)

        GTG = G.T.dot(G * w.reshape((-1, 1)))
        GTd = G.T.dot(w * d)
        GG = GTG + self.LTL * np.median(np.diagonal(GTG))

        try:
            c, lower = scipy.linalg.cho_factor(GG, lower=True)
            return scipy.linalg.cho_solve((c, lower), GTd)
        except scipy.linalg.LinAlgError:
            return scipy.linalg.lstsq(GG, GTd, cond=None, lapack_driver='gelsy')[0]
