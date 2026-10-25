Fixed consistent initial conditions of DAE models not being computed if
there was only a single output time point. Inconsistent initial values of
algebraic states were then reported as they were, or the simulation
failed.
