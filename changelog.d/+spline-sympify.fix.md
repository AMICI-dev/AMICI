Fixed two sympification gaps in the spline module (#2146, #3281):
spline node/value/derivative strings passed to the Python spline
construction API (`CubicHermiteSpline`, `UniformGrid`) went through a
bare `sympy.sympify()` with no `locals`, so identifiers could silently
collide with sympy builtins (e.g. a parameter named `I`) or fail to
compare equal to the same-named symbol used elsewhere in the model;
and `AbstractSpline.from_annotation` never substituted local symbols
into spline formulas parsed from `<amici:spline>` annotations, unlike
every other SBML math expression in the model.
