import numpy as np
import pytest

import cppe

from .utils import benchmark_settings, time_callable


_LARGE_POTFILE = "tests/potfiles/loprop_solvated_20.pot"
_NSITES = [32, 64, 128, 256]


def _subset_polsites(n_polsites):
    potentials = cppe.cpp.PotfileReader(_LARGE_POTFILE).read()
    polsites = cppe.get_polarizable_sites(potentials)
    if len(polsites) < n_polsites:
        raise RuntimeError(f"Need at least {n_polsites} polarizable sites.")
    return polsites[:n_polsites]


@pytest.mark.benchmark
@pytest.mark.parametrize("backend", ["cpp", "python"])
@pytest.mark.parametrize("n_polsites", _NSITES)
def test_bmatrix_apply_scaling(benchmark_recorder, backend, n_polsites):
    cppe.set_backend(backend)
    polsites = _subset_polsites(n_polsites)
    options = {
        "potfile": _LARGE_POTFILE,
        "summation_induced_fields": "direct",
        "tree_ncrit": 64,
        "tree_expansion_order": 5,
    }
    bmatrix = cppe.BMatrix(polsites, options)
    vec = np.linspace(1.0, 2.0, 3 * n_polsites)

    repeat, warmup = benchmark_settings()
    timings = time_callable(lambda: bmatrix.apply(vec), repeat=repeat, warmup=warmup)

    benchmark_recorder(
        benchmark="bmatrix_apply",
        backend=backend,
        dataset="loprop_solvated_20.pot",
        n_polsites=n_polsites,
        **timings,
    )

    cppe.set_backend("cpp")
