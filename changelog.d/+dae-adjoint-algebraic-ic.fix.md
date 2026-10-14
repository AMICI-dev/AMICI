Fixed adjoint sensitivity analysis failing for DAE models in which algebraic
states enter the differential equations. The algebraic components of the
adjoint state are now made consistent before the backward integration and
after each data point or event (#3343).
