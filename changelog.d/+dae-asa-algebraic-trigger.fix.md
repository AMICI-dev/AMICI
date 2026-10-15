Fixed adjoint sensitivities for DAE models with events whose trigger depends
on algebraic states, e.g., a `piecewise` in an algebraic rule with a condition
on an algebraic state. The contributions of the algebraic states to the
adjoint update at the event were lost (#3392).
