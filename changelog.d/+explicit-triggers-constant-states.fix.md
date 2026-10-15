Event triggers and `piecewise` conditions whose trigger times depend on
state variables with zero time derivative are now handled at their
precomputed trigger times, instead of by root-finding. This includes
parameters that are targets of event assignments, which AMICI converts to
state variables, e.g., the condition changes of PEtab v2 problems. The trigger
times are re-evaluated after each event. Previously, root-finding for such
triggers could produce wrong results: For example, for the PEtab v2 version
of the benchmark problem Isensee_JCB2018, a time-dependent stimulus was
switched off again immediately after it was switched on, because the solver
reported a spurious reverse crossing of the trigger function after the
switch (#3365).
