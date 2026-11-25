Fixed the forward sensitivities of event times for models with conservation
laws and root functions that depend on model expressions, such as PySB
models with `Piecewise` expressions of observables, or SBML models with
`piecewise` only in rate rules. The contribution via the model expressions
was counted twice. If a root function depended on a species eliminated by a
conservation law whose total depends on the parameters, the generated code did
not compile. Affected models need to be re-imported.
