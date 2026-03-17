import numpy as np

import cppe

from .cache import cache


def test_nuclear_fields_python_matches_cpp():
    mol = cache.molecule["pna"]
    potentials = cppe.cpp.PotfileReader("tests/potfiles/pna_6w.pot").read()

    cppe.set_backend("cpp")
    fields_cpp = cppe.NuclearFields(mol, potentials).compute()

    cppe.set_backend("python")
    fields_py = cppe.NuclearFields(mol, potentials).compute()

    cppe.set_backend("cpp")
    np.testing.assert_allclose(fields_py, fields_cpp, atol=1e-14)


def test_nuclear_gradient_python_matches_cpp():
    mol = cache.molecule["pna"]
    potentials = cppe.cpp.PotfileReader("tests/potfiles/pna_6w.pot").read()

    cppe.set_backend("cpp")
    grad_cpp = cppe.NuclearFields(mol, potentials).nuclear_gradient()

    cppe.set_backend("python")
    grad_py = cppe.NuclearFields(mol, potentials).nuclear_gradient()

    cppe.set_backend("cpp")
    np.testing.assert_allclose(grad_py, grad_cpp, atol=1e-13)


def test_nuclear_fields_numba_nopython_compiles():
    from cppe import fields

    mol = cache.molecule["pna"]
    potentials = cppe.cpp.PotfileReader("tests/potfiles/pna_6w.pot").read()

    cppe.set_backend("python")
    cppe.NuclearFields(mol, potentials).compute()
    cppe.NuclearFields(mol, potentials).nuclear_gradient()
    cppe.set_backend("cpp")

    assert fields._nuclear_fields_kernel.nopython_signatures
    assert fields._nuclear_gradient_kernel.nopython_signatures
