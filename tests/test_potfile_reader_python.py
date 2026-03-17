from pathlib import Path

import numpy as np
import pytest

import cppe


def _serialize_potential(potential):
    multipoles = []
    for multipole in potential.multipoles:
        multipoles.append((multipole.k, np.array(multipole.values, dtype=float)))

    polarizability = None
    if potential.is_polarizable:
        polarizability = np.array(potential.polarizability.values, dtype=float)

    return {
        "index": potential.index,
        "element": potential.element,
        "position": np.array(potential.position, dtype=float),
        "exclusions": list(potential.exclusions),
        "multipoles": multipoles,
        "polarizability": polarizability,
    }


@pytest.mark.parametrize(
    "filename",
    [
        "pna_6w.pot",
        "pna_6w_isopol.pot",
        "loprop_solvated_20.pot",
    ],
)
def test_python_potfile_reader_matches_cpp(filename):
    potfile = Path(__file__).parent / "potfiles" / filename

    cpp_potentials = cppe.PotfileReader(str(potfile)).read()
    cppe.set_backend("python")
    py_potentials = cppe.PotfileReader(str(potfile)).read()
    cppe.set_backend("cpp")

    assert len(cpp_potentials) == len(py_potentials)

    for cpp_potential, py_potential in zip(cpp_potentials, py_potentials):
        cpp_data = _serialize_potential(cpp_potential)
        py_data = _serialize_potential(py_potential)

        assert cpp_data["index"] == py_data["index"]
        assert cpp_data["element"] == py_data["element"]
        assert cpp_data["exclusions"] == py_data["exclusions"]
        np.testing.assert_allclose(
            cpp_data["position"], py_data["position"], atol=1e-14
        )

        assert len(cpp_data["multipoles"]) == len(py_data["multipoles"])
        for (cpp_k, cpp_values), (py_k, py_values) in zip(
            cpp_data["multipoles"], py_data["multipoles"]
        ):
            assert cpp_k == py_k
            np.testing.assert_allclose(cpp_values, py_values, atol=1e-14)

        if cpp_data["polarizability"] is None:
            assert py_data["polarizability"] is None
        else:
            np.testing.assert_allclose(
                cpp_data["polarizability"], py_data["polarizability"], atol=1e-14
            )


def test_python_potfile_reader_missing_file_raises():
    cppe.set_backend("python")
    with pytest.raises(RuntimeError, match="Potential file does not exist"):
        cppe.PotfileReader("does-not-exist.pot")
    cppe.set_backend("cpp")
