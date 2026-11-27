Fixed wrong event triggering for trigger or `piecewise` conditions within
the tolerances of their threshold, if the model has states on very different
scales, and endlessly retriggering events if an event moves a trigger that
is a difference of states, e.g., `x1 - x2 >= 0`, back by less than the
tolerances. The direction of such conditions and their uncertainty are now
computed from analytic derivatives of the trigger functions, which are
generated during model import. Models need to be re-imported.
Also fixed events with a trigger like `(s > c) && (u == 0)`, where `u` only
changes at events, triggering when an event assignment moves `s` below `c`.
