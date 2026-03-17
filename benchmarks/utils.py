import os
import statistics
import time


def benchmark_settings():
    repeat = int(os.environ.get("CPPE_BENCH_REPEAT", "3"))
    warmup = int(os.environ.get("CPPE_BENCH_WARMUP", "1"))
    return repeat, warmup


def time_callable(func, *, repeat, warmup):
    for _ in range(warmup):
        func()

    samples = []
    for _ in range(repeat):
        t0 = time.perf_counter()
        func()
        samples.append(time.perf_counter() - t0)

    mean_s = statistics.mean(samples)
    std_s = statistics.pstdev(samples) if len(samples) > 1 else 0.0
    return {
        "repeat": repeat,
        "warmup": warmup,
        "mean_s": mean_s,
        "std_s": std_s,
        "min_s": min(samples),
        "max_s": max(samples),
    }
