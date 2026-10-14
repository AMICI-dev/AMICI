Fixed uncompilable model code for noise formulas (sigmas) that reference
expressions, such as targets of SBML assignment rules. These are now replaced
by their definitions. Sigmas that, as a result, depend on quantities that are
not available for their evaluation (e.g., states that cannot be expressed in
terms of the respective observable) are now rejected at import with an
informative error message (#3384).
