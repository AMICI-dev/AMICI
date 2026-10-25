Fixed forward and adjoint sensitivities when an event assignment changes
a `piecewise` condition without triggering another event, e.g., when the
assignment moves a variable from exactly the threshold to below it (#3436).
The jump in the right-hand side was ignored. For DAEs, the algebraic states
are now also made consistent with the switched equations in this case.
