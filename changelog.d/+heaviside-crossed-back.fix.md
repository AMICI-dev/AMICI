Fixed wrong simulation results or failures if an event assignment reverts a
trigger or `piecewise` condition that changed at the same time, e.g.,
`at x < thr: x = 1.5 * x`, or moves it exactly to its threshold, e.g., for a
bouncing ball (#3425), and for events that move their own trigger time, e.g.,
`at time >= s: s = s + 1`, which only triggered once. Trigger and `piecewise`
conditions within the tolerances of their threshold are now resolved by the
direction in which they change, also at the start of a simulation or
pre-equilibration, and for conditions on algebraic states of DAEs. The
`ROOT_AFTER_REINIT` warning was replaced by `TRIGGER_WITHIN_TOLERANCE`. A
switch that flips back and forth at the same time, e.g.,
`x' = piecewise(-1, x > 0, 1)` at `x = 0`, and events that keep triggering
each other at the same time now fail instead of giving wrong results or not
terminating.
