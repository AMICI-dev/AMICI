AMICI no longer depends on Boost.Chrono. CPU times in `ReturnData`
(`cpu_time`, `cpu_time_b`, and others) are now always measured per thread,
using the operating system's thread CPU clock, also when Boost is not
installed.
