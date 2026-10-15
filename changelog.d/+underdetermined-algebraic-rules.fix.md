SBML import now fails with an informative error if the algebraic rules
do not determine all variables they reference (e.g., a non-constant
parameter in an algebraic rule). Previously, such models were imported,
and the simulation failed in `IDACalcIC`.
