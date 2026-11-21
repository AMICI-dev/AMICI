Fixed forward and adjoint sensitivities for DAE models with simultaneous
events if an event reads an algebraic state that depends on states changed by
an event applied before it. The sensitivities and time derivatives of the
algebraic states were not updated after event assignments (#3437). Affected models need to be re-imported.
