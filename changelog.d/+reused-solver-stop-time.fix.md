Fixed a failing post-equilibration-only simulation (`timepoints = [inf]`)
with a `Solver` that was previously used for an aborted simulation,
e.g., due to `max_steps`, if the aborted simulation ended before the new `t0`
([#3569](https://github.com/AMICI-dev/AMICI/issues/3569)).
