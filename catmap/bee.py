"""
Simplified version of ase.dft.bee

Input:
    32x1 array with non-selfconsitent energies.

Output:
    array with 2000 perturbations to the BEEF-vdW energy.

Warning!
    Only works with BEEF-vdW xc functional.

"""

import numpy as np
from ase.dft.pars_beefvdw import uiOmega as omega


class BEEFEnsemble:
    def __init__(self, size=2000, seed=0):
        self.size = size
        self.seed = seed

    def get_ensemble_perturbations(self, contribs):
        coefs = self.get_beefvdw_ensemble_coefs(self.size, self.seed)
        de = np.dot(coefs, contribs)
        return de

    def get_beefvdw_ensemble_coefs(self, size=2000, seed=0):
        """Perturbation coefficients of the BEEF-vdW ensemble"""
        assert np.shape(omega) == (31, 31)
        W, V, generator = self.eigendecomposition(omega, seed)
        RandV = generator.randn(31, size)
        coefs = []
        Q = np.dot(V, np.diag(np.sqrt(W)))
        for j in range(size):
            v = RandV[:, j]
            coefs.append(np.dot(Q, v)[:])
        ensemble_coefs = np.concatenate(coefs, axis=0).reshape(-1, 31)
        PBEc_ens = -ensemble_coefs[:, 30]
        return (np.vstack((ensemble_coefs.T, PBEc_ens))).T

    def eigendecomposition(self, omega, seed=0):
        u, s, v = np.linalg.svd(omega)  # unsafe: W, V = np.linalg.eig(omega)
        generator = np.random.RandomState(seed)
        return s, v.T, generator
