import numpy as np

import cppe


def test_multipole_fields_python_matches_cpp_direct():
    potentials = cppe.cpp.PotfileReader("tests/potfiles/pna_6w.pot").read()
    options = {
        "potfile": "tests/potfiles/pna_6w.pot",
        "summation_induced_fields": "direct",
    }

    cppe.set_backend("cpp")
    ref = cppe.MultipoleFields(potentials, options).compute()

    cppe.set_backend("python")
    got = cppe.MultipoleFields(potentials, options).compute()

    cppe.set_backend("cpp")
    np.testing.assert_allclose(got, ref, atol=1e-12)


def test_multipole_fields_numba_nopython_compiles():
    from cppe import multipole_fields

    potentials = cppe.cpp.PotfileReader("tests/potfiles/pna_6w.pot").read()
    options = {
        "potfile": "tests/potfiles/pna_6w.pot",
        "summation_induced_fields": "direct",
    }

    cppe.set_backend("python")
    cppe.MultipoleFields(potentials, options).compute()
    cppe.set_backend("cpp")

    assert multipole_fields._multipole_fields_direct_kernel.nopython_signatures
