The chain rule for the state derivatives of model expressions (`dwdx`),
used for the Jacobian and sensitivity equations, now only iterates over
chains of time-dependent expressions. Time-independent expressions, such as
initial assignments of parameters or expressions that select
experiment-specific values in PEtab v2 problems, no longer add recursion
levels to every Jacobian evaluation. Models need to be re-imported to
benefit from this.
