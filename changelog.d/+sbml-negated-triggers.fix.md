SBML import now supports negated comparisons in event triggers and
piecewise conditions, such as `!(x <= 3)` or `!(x <= 3 && y > 1)`.
Before, these raised an `SBMLException`. This requires sbmlmath>=TODO
([#3516](https://github.com/AMICI-dev/AMICI/issues/3516)).
