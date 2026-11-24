Fixed adjoint sensitivity analysis for DAE models if two discontinuities or
timepoints are closer than the solver's time resolution (about 1e-13 relative
to the time). It failed in `IDAGetQuadB`, or the gradient was wrong (#3411).
