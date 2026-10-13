Fixed the JAX backend attributing a root/discontinuity crossing to the
wrong Heaviside variable whenever a model had an event whose *solved*
trigger time still referenced a state (e.g. a trigger comparing `time`
directly to a state, `time >= x` — `sympy.solve` happily returns `[x]`
even though `x` is dynamic, not static). `DEModel._reorder_events()`
(which determines the physical event/Heaviside-array order) and the JAX
exporter's `iroot`/`eroot`/`ih`/`eh` classification (which determines
the root-detection order) used two different criteria for "does this
event have an explicit (precomputable) trigger time", which could
silently permute the two orderings relative to each other, causing
state updates to be gated by the wrong event's condition (#3286).
