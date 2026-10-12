Fixed DAE simulations failing or being inaccurate after events. After
event assignments, the solver was re-initialized with an invalid state
derivative, and algebraic states and derivatives were not made consistent
with the updated differential states. Furthermore, after the first event,
the post-event re-initialization step was needlessly repeated before every
subsequent output timepoint (#3267).
