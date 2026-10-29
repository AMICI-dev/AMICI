Handling of Discontinuities
===========================

This document provides guidance and rationale on the implementation of events
in AMICI. Events include any discontinuities in the right hand side of the
differential equation. There are three types of discontinuities:

- **Solution Jump Discontinuities** can be created by SBML events or delta
  functions in the right hand side.

- **Right-Hand-Side Jump Discontinuities** result in removable
  discontinuities in the solution and can be created by Piecewise,
  Heaviside functions and other logical operators in the right hand side.

- **Right-Hand-Side Removable Discontinuities** do not lead to
  discontinuities in the solution, but may lead to discontinuous higher
  order temporal derivatives and can be created by functions such as max or
  min in the right hand side.

Mathematical Considerations
---------------------------

A detailed mathematical description of the required sensitivity formulas is
provided in

* Fröhlich, F., Theis, F. J., Rädler, J. O., & Hasenauer, J. (2017).
  Parameter estimation for dynamical systems with discrete events and logical
  operations. Bioinformatics, 33(7), 1049-1056.
  doi:`10.1093/bioinformatics/btw764 <https://doi.org/10.1093/bioinformatics/btw764>`_.

Algorithmic Considerations
--------------------------

Solution Jump Discontinuities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SUNDIALS by itself does not support solution jump discontinuities. We
implement support by accessing private SUNDIALS API in
:cpp:func:`amici::Solver::reset_state`,
:cpp:func:`amici::Solver::reinit_post_process` and
:cpp:func:`amici::Solver::reinit_post_process_b`. These functions reset interval
variables to initial values to simulate a fresh integration start, but
keep/update the solution history, which is important for adjoint solutions.


Right-Hand-Side Jump Discontinuities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In principle these discontinuities do not need any special treatment, but
empirically, the solver may overstep or completely ignore the discontinuity,
leading to poor solution quality. This is particularly problematic when
step size is large and changes in step size, which can be caused by
parameter changes, inclusion of forward sensitivities or during backward
solves, may alter solutions in unexpected ways. Accordingly, finite
difference approximations, forward sensitivities as well as adjoint
sensitivities will yield poor derivative approximations.

To address these issues, we use the built-in rootfinding functionality in
SUNDIALS, which pauses the solver at the locations of discontinuities and
avoids overstepping or ignoring of discontinuities.

Root-finding is not required if the trigger time of a discontinuity can be
computed in advance, e.g., for `time >= t_on`. This applies if the trigger
time only depends on parameters, constants, expressions thereof, and state
variables with zero time derivative (e.g., parameters that are changed by
events, which are converted to state variables). For such triggers, AMICI
stops the integration exactly at the trigger time. The trigger times are
re-evaluated after each event.

Another difficulty comes with the evaluation of Heaviside functions. After
or during processing of discontinuities, Heaviside functions need to be
evaluated at the left and right hand limit of discontinuities.
This is challenging as the solver may slightly over- or understep the
discontinuity timepoint by a small epsilon and limits have to be correctly
computed in both forward and backward passes.

To address this issue, AMICI uses a vector of Heaviside helper variables `h`
that keeps track of the values of the Heaviside functions that have the
respective root function as argument. These will be automatically updated
during events and take either 0 or 1 values as appropriate pre/post event
limits.

The solver reports all zero crossings of the root functions. Crossings that
do not change the respective Heaviside helper variable are ignored, e.g., a
root that was just crossed and is reported once more right after the
discontinuity was handled, or a root function that was resolved by its
direction at a discontinuity (see below) and crosses zero only later.

Event assignments may change the root functions, which may trigger further
events or move the root functions of the current events back across zero.
Therefore, after each event assignment, the signs of the root functions are
compared to the Heaviside helper variables. A root function within its
uncertainty of zero is treated as follows: An event assignment that moves a
root function onto zero triggers the respective event, as its trigger is true
at zero. Afterwards, the direction in which the root function moves decides,
e.g., after an assignment exactly onto the threshold that reverses the
direction of the root function, or if it ended up slightly below zero due to
roundoff. The same holds for root functions that crossed zero at the current
timepoint, and for root functions at zero that were not changed by the event
assignments, whose direction may have changed nevertheless.
The uncertainty is estimated from the change of the
root function for a state perturbation within the integration tolerances
and within the time tolerance of the root-finding.
Likewise, at the start of a simulation or steady-state computation, the
Heaviside helper variables of root functions at zero are set by the direction
in which they move, unless the respective event triggers at the initial time.

In order to fully support SBML events and Piecewise functions, events only
trigger at zero crossings from negative to positive. For Heaviside
functions, two root functions are used to keep track of
Heaviside functions and two Heaviside function helper variables will be
created, where one corresponds to the value of `Heaviside(...)` and one
to the value of `1-Heaviside(...)`. To ensure that Heaviside functions are
correctly evaluated at the beginning of the simulation, Heaviside functions
are implement as unit steps that evaluate to `1` at `0`. The arguments of
Heaviside functions are normalized such that respective properties of
Piecewise functions are conserved for the first Heaviside function variable.
Accordingly, the value of of the second helper variable is incorrect when
simulation starts when the respective Heaviside function evaluates to zero
at initialization and should generally not be used.

States Around Discontinuities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At a discontinuity, several events may trigger at the same time and are
applied one after another, and their assignments may trigger further events.
The implementation distinguishes the following states, which all belong to
the same timepoint:

- ``*_at_trigger``: the state at which a discontinuity was reached, i.e.,
  when its events triggered, before the Heaviside helper variables were
  updated and before any event was applied. For events that are triggered
  by the assignments of other events, this is the state right after those
  assignments.
- ``*_pre`` / ``*_post``: the state right before / after the assignments of
  a single event were applied.
- ``*_eval``: the state at which the assignments of an event are evaluated.
  This is the ``*_at_trigger`` state for events that use values from the
  trigger time (SBML ``useValuesFromTriggerTime="true"``), and the ``*_pre``
  state otherwise.

Right-Hand-Side Removable Discontinuities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Removable discontinuities do not require any special treatment. Numerically,
this may be advantageous, but is currently not implemented.
