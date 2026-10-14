Fixed a crash when freeing an IDAS solver after an adjoint simulation of a DAE
model had failed before the first integration step (e.g., during
pre-equilibration by Newton's method). This is a patch of the bundled
SUNDIALS.
