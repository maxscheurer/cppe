import numpy as np
import pytest

import cppe

from .synthetic_systems import replicated_waterbox_system
from .utils import benchmark_settings, time_callable


_REPEATS = [
    (4, 4, 4),
    (6, 6, 6),
]


@pytest.mark.benchmark_large
@pytest.mark.parametrize("backend", ["cpp", "python"])
@pytest.mark.parametrize("repeats", _REPEATS)
def test_synthetic_bmatrix_apply(benchmark_recorder, backend, repeats):
    system = replicated_waterbox_system(repeats=repeats)
    options = {
        "potfile": "tests/potfiles/pna_6w.pot",
        "summation_induced_fields": "direct",
    }
    polsites = cppe.get_polarizable_sites(system["potentials"])
    vec = np.linspace(0.1, 0.7, 3 * len(polsites))

    cppe.set_backend(backend)
    bmat = cppe.BMatrix(polsites, options)
    repeat, warmup = benchmark_settings()
    timings = time_callable(lambda: bmat.apply(vec), repeat=repeat, warmup=warmup)

    benchmark_recorder(
        benchmark="synthetic_bmatrix_apply",
        backend=backend,
        dataset=system["label"],
        n_polsites=system["n_polsites"],
        **timings,
    )

    cppe.set_backend("cpp")


@pytest.mark.benchmark_large
@pytest.mark.parametrize("backend", ["cpp", "python"])
@pytest.mark.parametrize("repeats", _REPEATS)
def test_synthetic_nuclear_fields(benchmark_recorder, backend, repeats):
    system = replicated_waterbox_system(repeats=repeats)

    cppe.set_backend(backend)
    nf = cppe.NuclearFields(system["molecule"], system["potentials"])
    repeat, warmup = benchmark_settings()
    timings = time_callable(nf.compute, repeat=repeat, warmup=warmup)

    benchmark_recorder(
        benchmark="synthetic_nuclear_fields",
        backend=backend,
        dataset=system["label"],
        n_polsites=system["n_polsites"],
        **timings,
    )

    cppe.set_backend("cpp")
