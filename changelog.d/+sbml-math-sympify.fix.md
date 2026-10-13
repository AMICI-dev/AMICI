Some SBML-derived math (PEtab v1 initial assignments, assignment-rule-
based observables) was converted to L3 formula text and reparsed with
plain `sympy.sympify`, subject to sympy's own parsing grammar rather
than SBML's defined semantics. This math is now converted directly via
`sbmlmath`, like the rest of the SBML import pipeline; in rare cases,
some special constructs may now be interpreted differently (#2146,
#3281).
