from pathlib import Path

import cppe
import pytest

from .utils import benchmark_settings, time_callable


DATASETS = [
    "pna_6w.pot",
    "pna_6w_isopol.pot",
    "loprop_solvated_20.pot",
]


@pytest.mark.benchmark
@pytest.mark.parametrize("backend", ["cpp", "python"])
@pytest.mark.parametrize("dataset", DATASETS)
def test_potfile_reader(benchmark_recorder, backend, dataset):
    cppe.set_backend(backend)
    potfile = Path("tests") / "potfiles" / dataset

    def _run():
        return cppe.PotfileReader(str(potfile)).read()

    repeat, warmup = benchmark_settings()
    timings = time_callable(_run, repeat=repeat, warmup=warmup)
    n_polsites = len(cppe.get_polarizable_sites(_run()))

    benchmark_recorder(
        benchmark="potfile_read",
        backend=backend,
        dataset=dataset,
        n_polsites=n_polsites,
        **timings,
    )

    cppe.set_backend("cpp")
