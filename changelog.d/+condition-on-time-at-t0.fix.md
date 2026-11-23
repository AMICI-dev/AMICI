Fixed `piecewise` and trigger conditions on time that are at their threshold
at the initial time, e.g., `piecewise(2, time <= 0, 0)` simulated from
`t0 = 0`, which kept their value at `t0` for the whole simulation.
This includes trigger times that coincide with `t0` only up to roundoff,
e.g., `k * time >= 1` simulated from `t0 = 1 / k`
([#3515](https://github.com/AMICI-dev/AMICI/issues/3515)).
