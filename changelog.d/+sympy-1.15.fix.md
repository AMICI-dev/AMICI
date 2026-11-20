Fixed compatibility with SymPy 1.15: `And`/`Or`-combined event triggers
that SymPy now solves for `t` with undefined (`nan`) branches are no longer
treated as having explicit trigger times. With the JAX backend, this caused
time-triggered events of the affected models to be missed. With the SUNDIALS
backend, this resulted in spurious warnings about NaN trigger times.
