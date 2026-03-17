import csv
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import pytest

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("NUMBA_NUM_THREADS", "1")


def pytest_addoption(parser):
    parser.addoption(
        "--benchmark-output-dir",
        action="store",
        default="benchmarks/results",
        help="Directory for benchmark result artifacts.",
    )


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "benchmark: marks benchmark-style tests that collect timing metrics",
    )
    config.addinivalue_line(
        "markers",
        "benchmark_large: marks large synthetic benchmark tests",
    )


@pytest.fixture(scope="session", autouse=True)
def benchmark_thread_control():
    from numba import set_num_threads

    threads = int(os.environ.get("NUMBA_NUM_THREADS", "1"))
    set_num_threads(threads)


@pytest.fixture(scope="session")
def benchmark_recorder(request):
    records = []

    def _record(**kwargs):
        records.append(kwargs)

    yield _record

    output_dir = Path(request.config.getoption("benchmark_output_dir"))
    output_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    payload = {
        "created_utc": timestamp,
        "repeat": int(os.environ.get("CPPE_BENCH_REPEAT", "3")),
        "warmup": int(os.environ.get("CPPE_BENCH_WARMUP", "1")),
        "omp_num_threads": int(os.environ.get("OMP_NUM_THREADS", "1")),
        "numba_num_threads": int(os.environ.get("NUMBA_NUM_THREADS", "1")),
        "records": records,
    }

    latest_json = output_dir / "latest.json"
    latest_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    stamped_json = output_dir / f"benchmarks-{timestamp}.json"
    stamped_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    if records:
        fieldnames = sorted({k for rec in records for k in rec.keys()})
        latest_csv = output_dir / "latest.csv"
        stamped_csv = output_dir / f"benchmarks-{timestamp}.csv"
        for csv_path in (latest_csv, stamped_csv):
            with csv_path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(records)
