PEtab SciML: Neural network input expressions in hybridization tables are
now parsed with PEtab's math grammar, like observable formulas. Previously,
PEtab-specific functions such as `piecewise` were not recognized, and model
entities with IDs such as `lambda` could not be used in input expressions.
