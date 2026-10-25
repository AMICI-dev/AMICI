Fixed the results of simulations that failed before the main simulation
started, e.g., for unsupported settings (#3477).
`ReturnData` reported a log-likelihood, states, Jacobian and `xdot` evaluated
at zero states instead of NaN, and the CPU time, number of steps and final
time of the solver's previous simulation.
