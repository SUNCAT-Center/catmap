"""
Simplified version of ase.dft.bee. Only works with BEEF-vdW xc functional.

Input:
    32x1 array with non-selfconsistent energies.

Output:
    array with ensemble perturbations to the BEEF-vdW energy.
"""

import numpy as np
from ase.dft.pars_beefvdw import uiOmega as omega


class BEEFEnsemble:
    def __init__(self, size=2000, seed=0):
        """
        Parameters
        ----------
            size : int
                Length of the ensemble vector.
            seed : int
                Random seed. Must be identical for initial and final states,
                in order to calculate covariance. Therefore it is attached to
                the class.
        """

        self.size = size
        self.seed = seed

    def get_ellipse(self, de_x, de_y):
        """ Return width, height and angle of a BEEF covariance ellipse.

        Parameters
        ----------
            de_x : list
                Ensemble energy differences wrt. a fixed reference.
            de_y : list
                Ensemble energy differences wrt. a fixed reference.
        """
        cov = np.cov(de_x, de_y)
        eigval, eigvec = np.linalg.eig(cov)
        # Return the +/- 1 sigma confidence intervals.
        width, height = 2 * np.sqrt(eigval)
        angle = np.rad2deg(np.arccos(eigvec[0, 0]))
        return width, height, angle

    def get_ensemble_perturbations(self, contribs):
        """ Return BEEF ensemble.

        Parameters
        ----------
            contribs : list
                list of 32 non-selfconsistent energies from the
                BEEF-vdW exchange-correlation functional.
        """
        coefs = self.get_beefvdw_ensemble_coefs()
        ens = np.dot(coefs, contribs)
        return ens

    def get_beefvdw_ensemble_coefs(self):
        """Return perturbation coefficients of the BEEF-vdW ensemble.
        """
        assert np.shape(omega) == (31, 31)
        W, V, generator = self.eigendecomposition(omega)
        RandV = generator.randn(31, self.size)
        coefs = []
        Q = np.dot(V, np.diag(np.sqrt(W)))
        for j in range(self.size):
            v = RandV[:, j]
            coefs.append(np.dot(Q, v)[:])
        ensemble_coefs = np.concatenate(coefs, axis=0).reshape(-1, 31)
        PBEc_ens = -ensemble_coefs[:, 30]
        return (np.vstack((ensemble_coefs.T, PBEc_ens))).T

    def eigendecomposition(self, omega):
        u, s, v = np.linalg.svd(omega)
        generator = np.random.RandomState(self.seed)
        return s, v.T, generator
