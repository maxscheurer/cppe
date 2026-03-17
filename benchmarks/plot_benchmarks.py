#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def _parse_args():
    parser = argparse.ArgumentParser(
        description="Plot CPPE benchmark results and generate summary artifacts."
    )
    parser.add_argument(
        "--input",
        default="benchmarks/results/latest.csv",
        help="Path to benchmark CSV produced by benchmarks/conftest.py",
    )
    parser.add_argument(
        "--output-dir",
        default="benchmarks/plots",
        help="Directory for generated plots and summary tables.",
    )
    parser.add_argument(
        "--format",
        dest="formats",
        action="append",
        default=["png"],
        help="Output image format (repeatable), e.g. --format png --format svg",
    )
    parser.add_argument(
        "--dpi", type=int, default=160, help="Image DPI for raster plots."
    )
    parser.add_argument(
        "--no-log-y",
        action="store_true",
        help="Disable logarithmic y-axis on runtime plots.",
    )
    parser.add_argument(
        "--detailed",
        action="store_true",
        help="Also generate per-benchmark detailed plots.",
    )
    return parser.parse_args()


def _slug(text: str) -> str:
    return "".join(c if c.isalnum() else "_" for c in text).strip("_").lower()


def _save_fig(fig, out_base: Path, formats, dpi: int):
    for fmt in formats:
        fig.savefig(out_base.with_suffix(f".{fmt}"), dpi=dpi, bbox_inches="tight")


def _cleanup_detailed_outputs(out_dir: Path):
    for pattern in ("runtime_*", "speedup_*"):
        for path in out_dir.glob(pattern):
            if path.name.startswith("overview_"):
                continue
            if path.is_file():
                path.unlink()


def _aggregate(df: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        df.groupby(["benchmark", "dataset", "backend", "n_polsites"], as_index=False)
        .agg(
            mean_s=("mean_s", "mean"),
            std_s=("mean_s", "std"),
            min_s=("mean_s", "min"),
            max_s=("mean_s", "max"),
            samples=("mean_s", "count"),
        )
        .fillna({"std_s": 0.0})
    )
    return grouped


def _speedup_table(agg: pd.DataFrame) -> pd.DataFrame:
    pivot = agg.pivot_table(
        index=["benchmark", "dataset", "n_polsites"],
        columns="backend",
        values="mean_s",
    ).reset_index()
    if "cpp" in pivot.columns and "python" in pivot.columns:
        pivot["python_over_cpp"] = pivot["python"] / pivot["cpp"]
        pivot["cpp_over_python"] = pivot["cpp"] / pivot["python"]
    return pivot


def _plot_runtime_curves(
    agg: pd.DataFrame, out_dir: Path, formats, dpi: int, log_y: bool
):
    for (benchmark, dataset), sub in agg.groupby(["benchmark", "dataset"]):
        fig, ax = plt.subplots(figsize=(8.0, 5.0))
        for backend, bdf in sub.groupby("backend"):
            bdf = bdf.sort_values("n_polsites")
            ax.errorbar(
                bdf["n_polsites"],
                bdf["mean_s"],
                yerr=bdf["std_s"],
                marker="o",
                linewidth=1.8,
                capsize=3,
                label=backend,
            )

        ax.set_xscale("log", base=2)
        if log_y:
            ax.set_yscale("log")
        ax.set_xlabel("Number of polarizable sites (n_polsites)")
        ax.set_ylabel("Runtime [s]")
        ax.set_title(f"{benchmark} ({dataset})")
        ax.grid(alpha=0.3)
        ax.legend()

        out_base = out_dir / f"runtime_{_slug(benchmark)}_{_slug(dataset)}"
        _save_fig(fig, out_base, formats, dpi)
        plt.close(fig)


def _plot_speedups(speedup: pd.DataFrame, out_dir: Path, formats, dpi: int):
    needed = {"cpp", "python", "python_over_cpp"}
    if not needed.issubset(speedup.columns):
        return

    for (benchmark, dataset), sub in speedup.groupby(["benchmark", "dataset"]):
        sub = sub.sort_values("n_polsites")
        fig, ax = plt.subplots(figsize=(8.0, 4.5))
        ax.plot(
            sub["n_polsites"],
            sub["python_over_cpp"],
            marker="o",
            linewidth=1.8,
            color="#1f77b4",
        )
        ax.axhline(1.0, color="black", linestyle="--", linewidth=1.0)
        ax.set_xscale("log", base=2)
        ax.set_xlabel("Number of polarizable sites (n_polsites)")
        ax.set_ylabel("Python / C++ runtime")
        ax.set_title(f"Speedup ratio: {benchmark} ({dataset})")
        ax.grid(alpha=0.3)

        out_base = out_dir / f"speedup_{_slug(benchmark)}_{_slug(dataset)}"
        _save_fig(fig, out_base, formats, dpi)
        plt.close(fig)


def _plot_overview_runtime(
    agg: pd.DataFrame, out_dir: Path, formats, dpi: int, log_y: bool
):
    benchmarks = sorted(agg["benchmark"].unique())
    ncols = 2
    nrows = (len(benchmarks) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12, 4.4 * nrows))
    axes = axes.flatten() if hasattr(axes, "flatten") else [axes]

    for idx, benchmark in enumerate(benchmarks):
        ax = axes[idx]
        sub = agg[agg["benchmark"] == benchmark]
        for (backend, dataset), bdf in sub.groupby(["backend", "dataset"]):
            bdf = bdf.sort_values("n_polsites")
            label = (
                backend if sub["dataset"].nunique() == 1 else f"{backend} ({dataset})"
            )
            ax.plot(
                bdf["n_polsites"], bdf["mean_s"], marker="o", linewidth=1.7, label=label
            )

        ax.set_title(benchmark)
        ax.set_xscale("log", base=2)
        if log_y:
            ax.set_yscale("log")
        ax.set_xlabel("n_polsites")
        ax.set_ylabel("runtime [s]")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)

    for idx in range(len(benchmarks), len(axes)):
        fig.delaxes(axes[idx])

    fig.suptitle("CPPE Benchmark Runtime Overview", fontsize=14)
    _save_fig(fig, out_dir / "overview_runtime", formats, dpi)
    plt.close(fig)


def _plot_overview_speedup(speedup: pd.DataFrame, out_dir: Path, formats, dpi: int):
    needed = {"benchmark", "dataset", "n_polsites", "python_over_cpp"}
    if not needed.issubset(speedup.columns):
        return

    benchmarks = sorted(speedup["benchmark"].unique())
    ncols = 2
    nrows = (len(benchmarks) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12, 4.2 * nrows))
    axes = axes.flatten() if hasattr(axes, "flatten") else [axes]

    for idx, benchmark in enumerate(benchmarks):
        ax = axes[idx]
        sub = speedup[speedup["benchmark"] == benchmark]
        for dataset, bdf in sub.groupby("dataset"):
            bdf = bdf.sort_values("n_polsites")
            label = dataset if sub["dataset"].nunique() > 1 else "python/cpp"
            ax.plot(
                bdf["n_polsites"],
                bdf["python_over_cpp"],
                marker="o",
                linewidth=1.7,
                label=label,
            )

        ax.axhline(1.0, color="black", linestyle="--", linewidth=1.0)
        ax.set_title(benchmark)
        ax.set_xscale("log", base=2)
        ax.set_xlabel("n_polsites")
        ax.set_ylabel("python / cpp")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)

    for idx in range(len(benchmarks), len(axes)):
        fig.delaxes(axes[idx])

    fig.suptitle("CPPE Benchmark Speedup Overview", fontsize=14)
    _save_fig(fig, out_dir / "overview_speedup", formats, dpi)
    plt.close(fig)


def _write_report(out_dir: Path, agg: pd.DataFrame, speedup: pd.DataFrame):
    def _markdown_table(df: pd.DataFrame) -> str:
        headers = list(df.columns)
        sep = ["---"] * len(headers)
        rows = ["| " + " | ".join(headers) + " |", "| " + " | ".join(sep) + " |"]
        for _, row in df.iterrows():
            cells = [str(row[h]) for h in headers]
            rows.append("| " + " | ".join(cells) + " |")
        return "\n".join(rows)

    lines = [
        "# Benchmark Report",
        "",
        "## Artifacts",
        "- `summary_aggregated.csv`: aggregated runtime statistics",
        "- `summary_speedup.csv`: backend speedup ratios",
        "- `overview_runtime.*`: consolidated runtime overview",
        "- `overview_speedup.*`: consolidated speedup overview",
        "- `runtime_*.png|svg`: detailed runtime curves (with --detailed)",
        "- `speedup_*.png|svg`: detailed speedup curves (with --detailed)",
        "",
        "## Fastest mean runtime by benchmark/dataset",
        "",
    ]

    fastest = (
        agg.sort_values("mean_s")
        .groupby(["benchmark", "dataset"], as_index=False)
        .first()[["benchmark", "dataset", "backend", "n_polsites", "mean_s"]]
    )
    lines.append(_markdown_table(fastest))

    if {"cpp", "python", "python_over_cpp"}.issubset(speedup.columns):
        lines.extend(
            [
                "",
                "## Median Python/C++ ratio by benchmark",
                "",
            ]
        )
        ratio = (
            speedup.groupby(["benchmark", "dataset"], as_index=False)["python_over_cpp"]
            .median()
            .rename(columns={"python_over_cpp": "median_python_over_cpp"})
        )
        lines.append(_markdown_table(ratio))

    (out_dir / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    args = _parse_args()
    in_path = Path(args.input)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not args.detailed:
        _cleanup_detailed_outputs(out_dir)

    df = pd.read_csv(in_path)
    required = {"benchmark", "dataset", "backend", "n_polsites", "mean_s"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Input CSV is missing required columns: {sorted(missing)}")

    df["n_polsites"] = pd.to_numeric(df["n_polsites"])
    df["mean_s"] = pd.to_numeric(df["mean_s"])

    agg = _aggregate(df)
    speedup = _speedup_table(agg)

    agg.to_csv(out_dir / "summary_aggregated.csv", index=False)
    speedup.to_csv(out_dir / "summary_speedup.csv", index=False)

    _plot_overview_runtime(
        agg,
        out_dir,
        formats=args.formats,
        dpi=args.dpi,
        log_y=not args.no_log_y,
    )
    _plot_overview_speedup(speedup, out_dir, formats=args.formats, dpi=args.dpi)

    if args.detailed:
        _plot_runtime_curves(
            agg,
            out_dir,
            formats=args.formats,
            dpi=args.dpi,
            log_y=not args.no_log_y,
        )
        _plot_speedups(speedup, out_dir, formats=args.formats, dpi=args.dpi)
    _write_report(out_dir, agg, speedup)

    print(f"Wrote benchmark plots and summaries to: {out_dir}")


if __name__ == "__main__":
    main()
