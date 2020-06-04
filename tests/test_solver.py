import unittest
import os

import numpy as np

from cppe import CppeState, BMatrix

from cppe import get_polarizable_sites, InducedMoments
from cppe.tensors import T_recursive

from .test_functionality import print_callback
from .cache import cache

from scipy.sparse.linalg import LinearOperator


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
        options = {"potfile": self.potfile_path,
                   "induced_thresh": 1e-12}
        cppe_state = CppeState(options, mol, print_callback)
        cppe_state.calculate_static_energies_and_fields()

        static_fields = cppe_state.static_fields
        zeros = np.zeros_like(static_fields)
        cppe_state.update_induced_moments(zeros, False)

        potentials = cppe_state.potentials
        polsites = get_polarizable_sites(potentials)
        npolsites = cppe_state.get_polarizable_site_number()
        bmat = np.zeros((3 * npolsites, 3 * npolsites))

        for s1, pot1 in enumerate(polsites):
            pol = pot1.polarizability
            inv_alpha = triangle_to_mat(pol.values)
            bmat[block(s1, s1)] = np.linalg.inv(inv_alpha)
            for s2, pot2 in enumerate(polsites):
                if pot1.excludes_site(pot2.index) or s1 == s2:
                    continue
                diff = pot2.position - pot1.position
                T12 = T_recursive(2, diff)
                if s1 > s2:
                    bmat[block(s1, s2)] = -triangle_to_mat(T12)
                    bmat[block(s2, s1)] = -triangle_to_mat(T12)

        # Test the matrix apply
        bmatrix_cpp = BMatrix(polsites, options)
        ret = bmatrix_cpp.compute_apply(static_fields)
        ret_ref = bmat @ static_fields
        np.testing.assert_allclose(ret, ret_ref, atol=1e-12)

        ret1 = bmatrix_cpp.compute_apply_slice(static_fields, 0, 8)
        ret2 = bmatrix_cpp.compute_apply_slice(static_fields, 8, npolsites)
        ret_all = ret1 + ret2
        np.testing.assert_allclose(ret_all, ret_ref, atol=1e-12)

        # build the Bmatrix from C++
        A = LinearOperator(2 * (static_fields.size,),
                           matvec=bmatrix_cpp.compute_apply)
        Afull = A @ np.eye(static_fields.size)
        np.testing.assert_allclose(Afull, bmat, atol=1e-20)

        induced_moments_direct = np.linalg.inv(bmat) @ static_fields
        induced_moments_solver = cppe_state.get_induced_moments()
        np.testing.assert_almost_equal(induced_moments_solver,
                                       induced_moments_direct, decimal=9)

        ind_mom = InducedMoments(potentials, options)
        res_cg = ind_mom.compute_cg(static_fields)
        np.testing.assert_allclose(res_cg, induced_moments_direct, atol=1e-10)

        binv = bmatrix_cpp.direct_inverse()
        np.testing.assert_allclose(binv @ static_fields,
                                   induced_moments_direct, atol=1e-15)

    def test_solver_by_inversion_damped(self):
        mol = cache.molecule["pna"]
        options = {"potfile": self.potfile_path,
                   "damp_induced": True,
                   "induced_thresh": 1e-12}
        cppe_state = CppeState(options, mol, print_callback)
        cppe_state.calculate_static_energies_and_fields()

        static_fields = cppe_state.static_fields
        zeros = np.zeros_like(static_fields)
        cppe_state.update_induced_moments(zeros, False)

        potentials = cppe_state.potentials
        polsites = get_polarizable_sites(potentials)
        npolsites = cppe_state.get_polarizable_site_number()
        bmat = np.zeros((3 * npolsites, 3 * npolsites))

        for s1, pot1 in enumerate(polsites):
            pol = pot1.polarizability
            inv_alpha = triangle_to_mat(pol.values)
            bmat[block(s1, s1)] = np.linalg.inv(inv_alpha)
            for s2, pot2 in enumerate(polsites):
                if pot1.excludes_site(pot2.index) or s1 == s2:
                    continue
                diff = pot2.position - pot1.position
                alpha_i = pol.isotropic_value
                alpha_j = pot2.polarizability.isotropic_value
                v = 2.1304 / (alpha_i * alpha_j)**(1/6)
                T12 = T_recursive(2, diff, v)
                if s1 > s2:
                    bmat[block(s1, s2)] = -triangle_to_mat(T12)
                    bmat[block(s2, s1)] = -triangle_to_mat(T12)

        # Test the matrix apply
        bmatrix_cpp = BMatrix(polsites, options)
        ret = bmatrix_cpp.compute_apply(static_fields)
        ret_ref = bmat @ static_fields
        np.testing.assert_allclose(ret, ret_ref, atol=1e-12)

        ret1 = bmatrix_cpp.compute_apply_slice(static_fields, 0, 8)
        ret2 = bmatrix_cpp.compute_apply_slice(static_fields, 8, npolsites)
        ret_all = ret1 + ret2
        np.testing.assert_allclose(ret_all, ret_ref, atol=1e-12)

        # build the Bmatrix from C++
        A = LinearOperator(2 * (static_fields.size,),
                           matvec=bmatrix_cpp.compute_apply)
        Afull = A @ np.eye(static_fields.size)
        np.testing.assert_allclose(Afull, bmat, atol=1e-20)

        induced_moments_direct = np.linalg.inv(bmat) @ static_fields
        induced_moments_solver = cppe_state.get_induced_moments()
        np.testing.assert_almost_equal(induced_moments_solver,
                                       induced_moments_direct, decimal=9)

        ind_mom = InducedMoments(potentials, options)
        res_cg = ind_mom.compute_cg(static_fields)
        np.testing.assert_allclose(res_cg, induced_moments_direct, atol=1e-10)

        binv = bmatrix_cpp.direct_inverse()
        np.testing.assert_allclose(binv @ static_fields,
                                   induced_moments_direct, atol=1e-15)
