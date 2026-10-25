AMICI no longer depends on Boost.Chrono. On Linux and macOS, CPU times in
`ReturnData` (`cpu_time`, `cpu_time_b`, and others) are now always measured
per thread, using the operating system's thread CPU clock, also when Boost is
not installed. On Windows, whose per-thread CPU clock is too coarse for timing
single simulations, they are measured with `std::clock`, as without Boost
before.
