import pytest

import cppe

from tests.cache import cache

from .utils import benchmark_settings, time_callable


_LARGE_POTFILE = "tests/potfiles/loprop_solvated_20.pot"
_NSITES = [32, 64, 128, 256]


def _subset_potentials(n_polsites):
    potentials = cppe.cpp.PotfileReader(_LARGE_POTFILE).read()
    polsites = cppe.get_polarizable_sites(potentials)
    if len(polsites) < n_polsites:
        raise RuntimeError(f"Need at least {n_polsites} polarizable sites.")
    return polsites[:n_polsites]


@pytest.mark.benchmark
@pytest.mark.parametrize("backend", ["cpp", "python"])
@pytest.mark.parametrize("n_polsites", _NSITES)
def test_nuclear_fields_scaling(benchmark_recorder, backend, n_polsites):
    cppe.set_backend(backend)
    potentials = _subset_potentials(n_polsites)
    molecule = cache.molecule["nilered"]
    fields = cppe.NuclearFields(molecule, potentials)

    repeat, warmup = benchmark_settings()
    timings = time_callable(fields.compute, repeat=repeat, warmup=warmup)

    benchmark_recorder(
        benchmark="nuclear_fields",
        backend=backend,
        dataset="loprop_solvated_20.pot",
        n_polsites=n_polsites,
        **timings,
    )

    cppe.set_backend("cpp")
