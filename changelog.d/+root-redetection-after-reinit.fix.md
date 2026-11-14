Fixed `piecewise` switches that were reverted, and events that were
executed spuriously, right after a root if the solver stopped exactly at the
root, depending on the output timepoints (#3444).
