PEtab v1: Observable formulas are now parsed with PEtab's math grammar,
like noise formulas. Previously, they were re-parsed with `sympy.sympify`,
so PEtab-specific syntax such as `piecewise`, `&&`, `||`, or `!` caused
import errors, and `inf` was imported as an undefined symbol.
