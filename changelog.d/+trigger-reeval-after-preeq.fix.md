Fixed event triggers (including `piecewise` Heaviside functions) not being
re-evaluated from the re-initialized state at the start of the main
simulation (or pre-simulation) after pre-equilibration. If a trigger
depended on a state that was re-initialized by the simulation condition,
the trigger state from pre-equilibration was kept, and the subsequent
crossing was counted twice (#3316).
