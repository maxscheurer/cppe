import numpy as np

import cppe


def _polsites():
    potentials = cppe.cpp.PotfileReader("tests/potfiles/pna_6w.pot").read()
    return cppe.get_polarizable_sites(potentials)


def test_bmatrix_apply_python_matches_cpp_direct():
    polsites = _polsites()
    options = {
        "potfile": "tests/potfiles/pna_6w.pot",
        "summation_induced_fields": "direct",
    }
    vec = np.linspace(-0.2, 0.7, 3 * len(polsites))

    cppe.set_backend("cpp")
    ref = cppe.BMatrix(polsites, options).apply(vec)

    cppe.set_backend("python")
    got = cppe.BMatrix(polsites, options).apply(vec)
    cppe.set_backend("cpp")

    np.testing.assert_allclose(got, ref, atol=1e-12)


def test_bmatrix_dense_python_matches_cpp_direct():
    polsites = _polsites()
    options = {
        "potfile": "tests/potfiles/pna_6w.pot",
        "summation_induced_fields": "direct",
    }

    cppe.set_backend("cpp")
    ref = cppe.BMatrix(polsites, options).to_dense_matrix()

    cppe.set_backend("python")
    got = cppe.BMatrix(polsites, options).to_dense_matrix()
    cppe.set_backend("cpp")

    np.testing.assert_allclose(got, ref, atol=1e-12)


def test_bmatrix_python_numba_nopython_compiles():
    from cppe import solver

    polsites = _polsites()
    options = {
        "potfile": "tests/potfiles/pna_6w.pot",
        "summation_induced_fields": "direct",
    }
    vec = np.linspace(0.1, 1.1, 3 * len(polsites))

    cppe.set_backend("python")
    bmat = cppe.BMatrix(polsites, options)
    bmat.apply(vec)
    bmat.apply_diagonal(vec)
    bmat.apply_diagonal_inverse(vec)
    cppe.set_backend("cpp")

    assert solver._induced_fields_direct_kernel.nopython_signatures
    assert solver._apply_diagonal_kernel.nopython_signatures
    assert solver._apply_diagonal_inverse_kernel.nopython_signatures


def test_triangle_to_mat_matches_cpp_math_helper():
    from cppe import solver

    tri = np.array([1.1, -0.2, 0.7, 2.5, -0.9, 3.2])
    got = solver._triangle_to_mat(tri)
    ref = cppe.triangle_to_mat(tri.tolist())
    np.testing.assert_allclose(got, ref, atol=1e-14)
