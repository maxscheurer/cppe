import unittest
import numpy as np
from cppe import (Tk_tensor, xyz2idx,
                  factorial, prefactors, Tk_coefficients,
                  multipole_components)

from cppe import multipole_derivative
from . import symmetry_helper as symmetry

from polarizationsolver import fields as polfields
from polarizationsolver import tensors


prefs = [
    [-1.0],
    [-1.0, -1.0, -1.0],
    [-0.5, -1.0, -1.0, -0.5, -1.0, -0.5],
]


class TestMath(unittest.TestCase):

    def test_basic_math(self):
        # factorial
        rng = np.arange(0, 10, 1)
        for a in rng:
            assert factorial(a) == np.math.factorial(a)

    def test_prefactors(self):
        for k in range(3):
            assert prefactors(k) == prefs[k]

    def test_multipole_components(self):
        for k in range(6):
            comps = (k + 1) * (k + 2) / 2
            assert comps == multipole_components(k)

    def test_xyz2idx(self):
        for k in range(6):
            combinations = []
            for x in range(k, -1, -1):
                for y in range(k, -1, -1):
                    for z in range(k, -1, -1):
                        if (x + y + z) != k:
                            continue
                        combinations.append((x, y, z))
            for i, c in enumerate(combinations):
                assert xyz2idx(*c) == i

    def test_T_tensors(self):
        # tests the T tensors against auto-generated Python code
        for k in range(7):
            R = 2.0 * np.random.random(3) + 1
            ref = tensors.T[k](R)

            coeffs = Tk_coefficients(k)
            actual = Tk_tensor(k, R, coeffs)

            # gets the indices of non-redundant components
            # e.g., takes only xy from (xy, yx) and so on ...
            sym_indices = symmetry.get_symm_indices(k)
            np.testing.assert_almost_equal(
                actual, ref.take(sym_indices), decimal=10,
                err_msg="T tensors do not match. Order = {}".format(k)
            )

    def test_T_tensors_damped(self):
        # tests the T tensors against auto-generated Python code
        for k in range(4):
            R = 2.0 * np.random.random(3) + 1
            damp = 2.0
            a_i = 4.0
            a_j = 10.0
            a = 1/(a_i*a_j)**(1/6)*damp
            ref = tensors.T_damp_thole[k](R, a)

            coeffs = Tk_coefficients(k)
            # damped tensor: damping_factor, alpha_i, alpha_j
            actual = Tk_tensor(k, R, coeffs, damp, a_i, a_j)

            # gets the indices of non-redundant components
            # e.g., takes only xy from (xy, yx) and so on ...
            sym_indices = symmetry.get_symm_indices(k)
            np.testing.assert_almost_equal(
                actual, ref.take(sym_indices), decimal=10,
                err_msg="Damped T tensors do not match. Order = {}".format(k)
            )

    def test_multipole_derivative(self):
        Rab = 2.0 * np.random.random(3) + 1
        for k in range(6):
            ll = 1  # derivative
            n_comps_symm = multipole_components(k)
            M = np.random.random(n_comps_symm)

            # unfold the tensor
            M_full = symmetry.unfold_tensor(M, k)
            sym_indices = symmetry.get_symm_indices(k)
            np.testing.assert_allclose(M, M_full.take(sym_indices))

            coeffs = Tk_coefficients(k + ll + 1)
            res = multipole_derivative(k, ll, Rab, M, coeffs, 0.0, 0.0, 0.0)
            ref_field = polfields.field(Rab, k, M_full, ll)
            np.testing.assert_allclose(ref_field, res, atol=1e-14)

            if k <= 2:
                damp = 2.0
                a_i = 4.0
                a_j = 10.0
                a = 1/(a_i*a_j)**(1/6)*damp
                res = multipole_derivative(k, ll, Rab, M, coeffs,
                                           damp, a_i, a_j)
                ref_field = polfields.thole_exp_field(Rab, k, M_full, ll, a)
                np.testing.assert_allclose(ref_field, res, atol=1e-14)
