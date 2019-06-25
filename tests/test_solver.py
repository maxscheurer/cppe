import unittest
import os
import h5py

import numpy as np

from cppe import Atom, Molecule
from cppe import PeOptions, CppeState

from cppe import (Tk_tensor,
                  Tk_coefficients,
                  get_polarizable_sites)

from .test_functionality import print_callback
from .cache import cache


def triangle_to_mat(y):
    mat = np.empty((3, 3))
    mat[0, 0] = y[0]
    mat[0, 1] = mat[1, 0] = y[1]
    mat[0, 2] = mat[2, 0] = y[2]
    mat[1, 1] = y[3]
    mat[1, 2] = mat[2, 1] = y[4]
    mat[2, 2] = y[5]
    return mat


def block(s1, s2):
    return (slice(3 * s1, 3 * s1 + 3), slice(3 * s2, 3 * s2 + 3))


class TestSolver(unittest.TestCase):
    dirname = os.path.dirname(__file__)
    potfile_path = "{}/potfiles/pna_6w.pot".format(dirname)

    def test_solver_by_inversion(self):
        mol = cache.molecule["pna"]
        options = PeOptions()
        options.potfile = self.potfile_path
        options.induced_thresh = 12
        options.diis_start_norm = 10
        cppe_state = CppeState(options, mol, print_callback)
        cppe_state.calculate_static_energies_and_fields()

        static_fields = np.array(cppe_state.get_static_fields())
        zeros = np.zeros_like(static_fields)
        cppe_state.update_induced_moments(zeros, False)

        potentials = cppe_state.get_potentials()
        polsites = get_polarizable_sites(potentials)
        npolsites = cppe_state.get_polarizable_site_number()
        bmat = np.zeros((3 * npolsites, 3 * npolsites))

        coeffs = Tk_coefficients(5)
        for s1, pot1 in enumerate(polsites):
            for pol in pot1.polarizabilities:
                inv_alpha = triangle_to_mat(pol.values)
                bmat[block(s1, s1)] = np.linalg.inv(inv_alpha)
            for s2, pot2 in enumerate(polsites):
                if pot1.excludes_site(pot2.index) or s1 == s2:
                    continue
                diff = pot2.position - pot1.position
                T12 = Tk_tensor(2, diff, coeffs)
                if s1 > s2:
                    bmat[block(s1, s2)] = -triangle_to_mat(T12)
                    bmat[block(s2, s1)] = -triangle_to_mat(T12)
        induced_moments_direct = np.linalg.inv(bmat) @ static_fields
        induced_moments_solver = cppe_state.get_induced_moments()
        np.testing.assert_almost_equal(induced_moments_solver,
                                       induced_moments_direct, decimal=9)
