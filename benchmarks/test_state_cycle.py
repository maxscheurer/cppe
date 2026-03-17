import numpy as np
import pytest

import cppe

from tests.cache import cache

from .utils import benchmark_settings, time_callable


@pytest.mark.benchmark
@pytest.mark.parametrize("backend", ["cpp", "python"])
@pytest.mark.parametrize("dataset", ["pna_6w.pot", "pna_6w_isopol.pot"])
def test_state_cycle(benchmark_recorder, backend, dataset):
    cppe.set_backend(backend)
    options = {"potfile": f"tests/potfiles/{dataset}", "induced_thresh": 1e-10}
    molecule = cache.molecule["pna"]

    def _run():
        state = cppe.CppeState(options, molecule, lambda _: None)
        state.calculate_static_energies_and_fields()
        zeros = np.zeros_like(state.static_fields)
        state.update_induced_moments(zeros, False)
        return state

    repeat, warmup = benchmark_settings()
    timings = time_callable(_run, repeat=repeat, warmup=warmup)
    n_polsites = _run().get_polarizable_site_number()

    benchmark_recorder(
        benchmark="state_cycle",
        backend=backend,
        dataset=dataset,
        n_polsites=n_polsites,
        **timings,
    )

    cppe.set_backend("cpp")
