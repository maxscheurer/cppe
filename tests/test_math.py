import unittest
import itertools
import numpy as np
from cppe import smat_vec, Tk_tensor, xyz2idx, factorial, prefactors, Tk_coefficients

from tensors import tensors
from tensors import symmetry


  # m.def("multipole_derivative", &multipole_derivative);
  # m.def("T", &T);
  # m.def("Tk_tensor", &Tk_tensor);
  # m.def("Tk_coefficients", &Tk_coefficients);
  # m.def("xyz2idx", &xyz2idx);
  # m.def("factorial", &factorial);
  # m.def("make_df", &make_df);
  # m.def("trinom", &trinom);
  # m.def("symmetry_factors", &symmetry_factors);
  # m.def("multipole_components", &multipole_components);

prefs = [
    [-1.0],
    [-1.0, -1.0, -1.0],
    [-0.5, -1.0, -1.0, -0.5, -1.0, -0.5],
]


class TestMath(unittest.TestCase):

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

    def test_T_tensors(self):
        ref_T = tensors.T

        for k in range(6):
            R = np.random.random(3)
            ref = ref_T[k](*R)

            coeffs = Tk_coefficients(k)
            actual = Tk_tensor(k, R, coeffs)

            sym_indices = symmetry.get_symm_indices(k)
            np.testing.assert_almost_equal(actual,
                                           ref.take(sym_indices), decimal=11)
