Fixed forward sensitivity analysis of DAE models failing in `IDACalcIC`
for some models with the default sensitivity tolerances. The initial
derivatives of the differential states and their sensitivities are now
computed from the model equations instead of being initialized to zero.
