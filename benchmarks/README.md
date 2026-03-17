# CPPE Benchmarks

This directory contains pytest-based benchmark suites that compare CPPE backends
and track scaling with the number of polarizable sites.

## Run

From a configured development environment:

```bash
python -m pytest benchmarks -m benchmark
```

Single-thread parity (recommended for fair comparisons):

```bash
OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 python -m pytest benchmarks -m benchmark
```

Optional runtime controls:

- `CPPE_BENCH_REPEAT` (default: `3`)
- `CPPE_BENCH_WARMUP` (default: `1`)
- `--benchmark-output-dir` (default: `benchmarks/results`)

Example:

```bash
CPPE_BENCH_REPEAT=5 CPPE_BENCH_WARMUP=2 \
python -m pytest benchmarks -m benchmark --benchmark-output-dir benchmarks/results
```

## Outputs

Each run writes both JSON and CSV artifacts:

- `benchmarks/results/latest.json`
- `benchmarks/results/latest.csv`
- timestamped copies for history

## Plotting

Generate runtime and speedup plots from benchmark CSV data:

```bash
python benchmarks/plot_benchmarks.py \
  --input benchmarks/results/latest.csv \
  --output-dir benchmarks/plots \
  --format png
```

By default, this writes **compact overview plots** only.

Optional:

- add additional output formats: `--format svg`
- disable log y-axis: `--no-log-y`
- generate per-benchmark detailed plots: `--detailed`

The script writes:

- `benchmarks/plots/summary_aggregated.csv`
- `benchmarks/plots/summary_speedup.csv`
- `benchmarks/plots/overview_runtime.*`
- `benchmarks/plots/overview_speedup.*`
- `benchmarks/plots/runtime_*.png` and `speedup_*.png` (only with `--detailed`)
- `benchmarks/plots/REPORT.md`

## Large synthetic systems

Large synthetic waterbox-like benchmarks are provided separately and are not run in
the default benchmark marker selection.

```bash
OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 python -m pytest benchmarks -m benchmark_large
```
