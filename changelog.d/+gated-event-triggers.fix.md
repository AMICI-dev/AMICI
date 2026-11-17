Fixed SBML event triggers with conditions on quantities that only change at events
(parameters, or species and parameters that are only changed by events).
Combined with another condition, e.g., `(x > c) && (u == 0)`,
such events were triggered at the end of the next integration step,
often the next output time, instead of when `x` crossed `c` (#3586).
Strict comparisons and `!=` on such quantities at their threshold were treated as true,
e.g., `(p > 0) && (time > 10)` with `p = 0` triggered at `t = 10` (#3508).
These conditions are now evaluated exactly instead of being part of the root function.
Also, time-dependent triggers containing functions that don't depend on time,
e.g., `time >= exp(p)`, were treated as periodic.
Their trigger times are now determined symbolically where possible,
instead of by root-finding.
