import unittest
import numpy as np
from cppe import (smat_vec, Tk_tensor, xyz2idx,
                  factorial, prefactors, Tk_coefficients,
                  multipole_components)

from .tensors import tensors, symmetry

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

    def test_smat_vec(self):
        x = np.random.random(3)
        y = np.random.random(6)
        z = smat_vec(y, x, 1.0)
        mat = np.empty((3, 3))
        mat[0, 0] = y[0]
        mat[0, 1] = mat[1, 0] = y[1]
        mat[0, 2] = mat[2, 0] = y[2]
        mat[1, 1] = y[3]
        mat[1, 2] = mat[2, 1] = y[4]
        mat[2, 2] = y[5]
        np.testing.assert_almost_equal(mat @ x, z, decimal=14)

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
        ref_T = tensors.T

        for k in range(6):
            R = np.random.random(3)
            ref = ref_T[k](*R)

            coeffs = Tk_coefficients(k)
            actual = Tk_tensor(k, R, coeffs)

            # gets the indices of non-redundant components
            # e.g., takes only xy from (xy, yx) and so on ...
            sym_indices = symmetry.get_symm_indices(k)
            np.testing.assert_almost_equal(actual,
                                           ref.take(sym_indices), decimal=10)

    def test_T_tensors_damped(self):
        # tests the T tensors against auto-generated Python code
        ref_T = tensors.T_damp_thole

        for k in range(4):
            R = np.random.random(3)
            damp = 2.0
            a_i = 4.0
            a_j = 10.0
            a = 1/(a_i*a_j)**(1/6)*damp
            ref = ref_T[k](*R, a)

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
