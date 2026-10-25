`amici::AbstractModel`'s virtual methods `fdJydsigma`, `fdJzdsigma`,
`fdJrzdsigma`, and `fdJrzdz` have been renamed to `fdJydsigmay`,
`fdJzdsigmaz`, `fdJrzdsigmaz`, and `fdJrzdrz` respectively, to correctly
reflect the symbol each derivative is actually taken with respect to
(`sigmay`, `sigmaz`, `rz`) rather than a generic, imprecise `sigma`/`z`
inherited from legacy MATLAB-era code generation. This only affects users
who subclass `AbstractModel`/`amici::Model` directly in C++; generated
model code from the standard SBML/PySB Python code generator is
unaffected other than being regenerated with the new names.
