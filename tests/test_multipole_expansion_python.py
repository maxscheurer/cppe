import numpy as np

import cppe

from .cache import cache


def test_multipole_interaction_energy_python_matches_cpp():
    mol = cache.molecule["pna"]
    potentials = cppe.cpp.PotfileReader("tests/potfiles/pna_6w.pot").read()

    cppe.set_backend("cpp")
    e_cpp = cppe.MultipoleExpansion(mol, potentials).interaction_energy()

    cppe.set_backend("python")
    e_py = cppe.MultipoleExpansion(mol, potentials).interaction_energy()

    cppe.set_backend("cpp")
    np.testing.assert_allclose(e_py, e_cpp, atol=1e-12)


def test_multipole_nuclear_gradient_python_matches_cpp():
    mol = cache.molecule["pna"]
    potentials = cppe.cpp.PotfileReader("tests/potfiles/pna_6w.pot").read()

    cppe.set_backend("cpp")
    grad_cpp = cppe.MultipoleExpansion(mol, potentials).nuclear_gradient()

    cppe.set_backend("python")
    grad_py = cppe.MultipoleExpansion(mol, potentials).nuclear_gradient()

    cppe.set_backend("cpp")
    np.testing.assert_allclose(grad_py, grad_cpp, atol=1e-12)


def test_multipole_interaction_numba_nopython_compiles():
    from cppe import multipole

    mol = cache.molecule["pna"]
    potentials = cppe.cpp.PotfileReader("tests/potfiles/pna_6w.pot").read()

    cppe.set_backend("python")
    cppe.MultipoleExpansion(mol, potentials).interaction_energy()
    cppe.set_backend("cpp")

    assert multipole._interaction_energy_kernel.nopython_signatures
