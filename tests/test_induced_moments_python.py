import numpy as np

import cppe

from .cache import cache


def test_induced_moments_python_matches_cpp_direct():
    mol = cache.molecule["pna"]
    options = {"potfile": "tests/potfiles/pna_6w.pot", "induced_thresh": 1e-12}

    cppe.set_backend("cpp")
    state = cppe.CppeState(options, mol, lambda _: None)
    state.calculate_static_energies_and_fields()
    rhs = state.static_fields
    potentials = state.potentials

    cppe.set_backend("cpp")
    ref = cppe.InducedMoments(potentials, options).compute(rhs, True)

    cppe.set_backend("python")
    got = cppe.InducedMoments(potentials, options).compute(rhs, True)

    cppe.set_backend("cpp")
    np.testing.assert_allclose(got, ref, atol=1e-10)
