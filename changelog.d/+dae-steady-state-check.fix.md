Fixed the steady-state computation for DAE models. The convergence checks
used the DAE residual, which vanishes along any consistent trajectory, so
that the integration-based steady-state search stopped right at the initial
state if no forward sensitivities were computed. The residual is now
evaluated at a zero state derivative (#3345).
