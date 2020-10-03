# Python 3.6.1

import numpy as np

from .base import AnalysisBase


class CurvatureAnalysis(AnalysisBase):
    def __init__(self, coords, shape, **kwargs):
        """
        :param coords: coords discreetly describing a surface
        :param shape:
        """
        super().__init__(coords, **kwargs)
        self.shape = shape

    def curvature(self):
        """
        :return: principal curvature 1, principal curvature 2, gaussian curvature, mean curvature
        """
        ll = len(self.coords)
        x = self.coords[:, 0].reshape(self.shape)
        y = self.coords[:, 1].reshape(self.shape)
        z = self.coords[:, 2].reshape(self.shape)

        # First derivatives
        x_v, x_u = np.gradient(x)
        y_v, y_u = np.gradient(y)
        z_v, z_u = np.gradient(z)

        # Second Derivatives
        x_vv, x_uv = np.gradient(x_v)
        x_uv, x_uu = np.gradient(x_u)
        y_vv, y_uv = np.gradient(y_v)
        y_uv, y_uu = np.gradient(y_u)
        z_vv, z_uv = np.gradient(z_v)
        z_uv, z_uu = np.gradient(z_u)

        # 2D to 1D for convenience
        x_v = np.reshape(x_v, ll)
        x_u = np.reshape(x_u, ll)
        y_v = np.reshape(y_v, ll)
        y_u = np.reshape(y_u, ll)
        z_v = np.reshape(z_v, ll)
        z_u = np.reshape(z_u, ll)
        x_vv = np.reshape(x_vv, ll)
        x_uv = np.reshape(x_uv, ll)
        x_uu = np.reshape(x_uu, ll)
        y_vv = np.reshape(y_vv, ll)
        y_uv = np.reshape(y_uv, ll)
        y_uu = np.reshape(y_uu, ll)
        z_vv = np.reshape(z_vv, ll)
        z_uv = np.reshape(z_uv, ll)
        z_uu = np.reshape(z_uu, ll)

        # Concatenate 1st and 2nd derivatives
        d1_v = np.c_[x_v, y_v, z_v]
        d1_u = np.c_[x_u, y_u, z_u]
        d2_vv = np.c_[x_vv, y_vv, z_vv]
        d2_uv = np.c_[x_uv, y_uv, z_uv]
        d2_uu = np.c_[x_uu, y_uu, z_uu]

        # 1st Fundamental Coefficients (E, F, G)
        e = np.sum(d1_v * d1_v, axis=1)
        f = np.sum(d1_u * d1_v, axis=1)
        g = np.sum(d1_u * d1_u, axis=1)

        # 2nd Fundamental Coefficients (L, M, N)
        cross = np.cross(d1_u, d1_v)
        norm = np.sqrt(np.sum(cross * cross, axis=1))
        unit = cross / np.c_[norm, norm, norm]
        l = np.sum(d2_vv * unit, axis=1)
        m = np.sum(d2_uv * unit, axis=1)
        n = np.sum(d2_uu * unit, axis=1)

        # Gaussian Curvature
        # K = det(second fundamental) / det(first fundamental)
        k = (l * n - m**2) / (e * g - f**2)

        # Mean Curvature
        # H = Tr[(second fundamental)(first fundamental inverse)] / 2
        h2 = np.trace(np.array([[l*g-f*m, -l*f+e*m], [m*g-n*f, -m*f+n*e]]) / (e*g-f**2))

        # Principle Curvatures
        diff = np.sqrt(abs(h2 ** 2 - 4 * k))
        kappa1 = h2 + diff
        kappa2 = h2 - diff

        return np.vstack([kappa1, kappa2, k, h2 / 2]).T
