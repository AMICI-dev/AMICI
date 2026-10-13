PEtab v1 import now parses observable and noise formulas using PEtab
v2's math grammar instead of a plain, ill-defined `sympy.sympify` call.
PEtab v1 never formally specified a grammar for these expressions; in
rare cases, some formulas may now be interpreted differently (#3281).
