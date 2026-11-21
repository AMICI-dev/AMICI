Added *static expressions* (`ws`) as a new model quantity
(`amici._symbolic.de_model_components.StaticExpression`): derived
quantities that are evaluated once at the initial time point from free
parameters, fixed parameters, the initial time, and other static
expressions. Unlike regular expressions (`w`), they can be referenced from
initial conditions, sigma formulas, and event assignments, and their
parameter dependence is accounted for in all sensitivities. Entries that
depend on fixed parameters are re-evaluated upon fixed-parameter-dependent
reinitialization after pre-equilibration; all other entries keep their
initial value. Their values are appended to `ReturnData::w` (shape
`nt` x (`nw` + `nws`)) and their IDs to the expression IDs. This is the
infrastructure for a simpler mapping of SBML initial assignments (#3409);
no importer makes use of it yet.

As a consequence, all generated model functions that take the fixed
parameters `k` now also take `ws`, and the corresponding
`amici::AbstractModel` virtual methods gained a `realtype const* ws`
argument. This only affects users who subclass `AbstractModel`/
`amici::Model` directly in C++; generated model code is unaffected other
than being regenerated with the new signatures.
