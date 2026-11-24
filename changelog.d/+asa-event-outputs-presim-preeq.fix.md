Fixed adjoint sensitivities for models with event outputs if events trigger
during pre-equilibration, pre-simulation, or post-equilibration: terms of
the main simulation's event likelihood were wrongly applied to these
periods ([#3413](https://github.com/AMICI-dev/AMICI/issues/3413)).
