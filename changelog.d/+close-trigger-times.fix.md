Fixed simulation failures with `CVODE routine reInitPostProcess failed with
error code -27` for events whose trigger times AMICI computes in advance
(e.g., `time >= p1 + p2`) if such a trigger time is within floating-point
roundoff of another one, of an output timepoint, or of the initial time,
e.g., for `time >= p1 + p2` and `time >= p3` with `p1 = 0.1`, `p2 = 0.2`,
and `p3 = 0.3`. Such events are now handled at the time at which the solver
stopped. An event that triggers within roundoff after an output timepoint is
therefore reflected in the output at that timepoint, as for an event that
triggers exactly at that timepoint. This does not apply to the last finite
output timepoint, after which such events are left to post-equilibration
(#3393).
