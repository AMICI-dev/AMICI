Fixed a check point being created for adjoint sensitivity analysis of DAE
models if the forward integration took exactly `maxsteps` steps. As for ODE
models (#2645), the check point interval is now `maxsteps + 1`.
