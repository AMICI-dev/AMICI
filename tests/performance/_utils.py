"""Shared utilities for JAX performance tests."""

import time

import jax

N_TIMING_RUNS = 5
N_WARMUP_RUNS = 2


def measure_exec_time(fn, *args, **kwargs) -> tuple[float, float]:
    """Return (t_first_s, t_exec_s).

    ``t_first_s`` is the wall time of the very first call (may include JIT
    compilation if this is the first call with these argument shapes).
    ``t_exec_s`` is the minimum over N_TIMING_RUNS subsequent calls.
    """
    t0 = time.perf_counter()
    fn(*args, **kwargs)
    jax.effects_barrier()
    t_first = time.perf_counter() - t0

    for _ in range(N_WARMUP_RUNS - 1):
        fn(*args, **kwargs)
        jax.effects_barrier()

    times = []
    for _ in range(N_TIMING_RUNS):
        t0 = time.perf_counter()
        fn(*args, **kwargs)
        jax.effects_barrier()
        times.append(time.perf_counter() - t0)

    return t_first, min(times)
