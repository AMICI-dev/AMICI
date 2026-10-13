Fixed an out-of-bounds read, typically causing a crash, when computing
forward sensitivities with `InternalSensitivityMethod.staggered` for DAE
models. The nonlinear solver for the sensitivities was created for one
vector too many.
