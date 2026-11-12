Fixed the SBML import of models with `rateOf` in initial assignments, where
the rate depends on species. Before, the generated `x0` code referenced
undeclared state variables and the model failed to compile
([#3504](https://github.com/AMICI-dev/AMICI/issues/3504)).
