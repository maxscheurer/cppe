import numpy as np

import cppe

from .cache import cache


def test_state_static_fields_python_matches_cpp():
    mol = cache.molecule["pna"]
    options = {"potfile": "tests/potfiles/pna_6w.pot"}

    cppe.set_backend("cpp")
    state_cpp = cppe.CppeState(options, mol, lambda _: None)
    state_cpp.calculate_static_energies_and_fields()

    cppe.set_backend("python")
    state_py = cppe.CppeState(options, mol, lambda _: None)
    state_py.calculate_static_energies_and_fields()

    cppe.set_backend("cpp")
    np.testing.assert_allclose(
        state_py.nuclear_fields, state_cpp.nuclear_fields, atol=1e-13
    )
    np.testing.assert_allclose(
        state_py.multipole_fields, state_cpp.multipole_fields, atol=1e-13
    )
    np.testing.assert_allclose(
        state_py.static_fields, state_cpp.static_fields, atol=1e-13
    )
    np.testing.assert_allclose(
        state_py.energies["Electrostatic"]["Nuclear"],
        state_cpp.energies["Electrostatic"]["Nuclear"],
        atol=1e-12,
    )


def test_state_gradient_methods_python_matches_cpp():
    mol = cache.molecule["pna"]
    options = {"potfile": "tests/potfiles/pna_6w.pot"}

    cppe.set_backend("cpp")
    state_cpp = cppe.CppeState(options, mol, lambda _: None)
    grad_e_cpp = state_cpp.nuclear_interaction_energy_gradient()
    grad_f_cpp = state_cpp.nuclear_field_gradient()

    cppe.set_backend("python")
    state_py = cppe.CppeState(options, mol, lambda _: None)
    grad_e_py = state_py.nuclear_interaction_energy_gradient()
    grad_f_py = state_py.nuclear_field_gradient()

    cppe.set_backend("cpp")
    np.testing.assert_allclose(grad_e_py, grad_e_cpp, atol=1e-12)
    np.testing.assert_allclose(grad_f_py, grad_f_cpp, atol=1e-12)


def test_state_update_induced_moments_python_matches_cpp():
    mol = cache.molecule["pna"]
    options = {"potfile": "tests/potfiles/pna_6w.pot", "induced_thresh": 1e-12}

    cppe.set_backend("cpp")
    state_cpp = cppe.CppeState(options, mol, lambda _: None)
    state_cpp.calculate_static_energies_and_fields()
    zeros_cpp = np.zeros_like(state_cpp.static_fields)
    state_cpp.update_induced_moments(zeros_cpp, False)

    cppe.set_backend("python")
    state_py = cppe.CppeState(options, mol, lambda _: None)
    state_py.calculate_static_energies_and_fields()
    zeros_py = np.zeros_like(state_py.static_fields)
    state_py.update_induced_moments(zeros_py, False)

    cppe.set_backend("cpp")
    np.testing.assert_allclose(
        state_py.get_induced_moments(), state_cpp.get_induced_moments(), atol=1e-10
    )
    np.testing.assert_allclose(
        state_py.energies["Polarization"]["Electronic"],
        state_cpp.energies["Polarization"]["Electronic"],
        atol=1e-12,
    )
    np.testing.assert_allclose(
        state_py.energies["Polarization"]["Nuclear"],
        state_cpp.energies["Polarization"]["Nuclear"],
        atol=1e-12,
    )
    np.testing.assert_allclose(
        state_py.energies["Polarization"]["Multipoles"],
        state_cpp.energies["Polarization"]["Multipoles"],
        atol=1e-12,
    )
