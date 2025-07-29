.. _python_caveats:

=======
Caveats
=======

This section lists some caveats that are important to be aware of when using
AMICI.

See also the `list of known bugs <https://github.com/AMICI-dev/AMICI/issues?q=is%3Aissue%20state%3Aopen%20label%3Abug>`__.

Sensitivity analysis
====================

AMICI provides high-performance sensitivity analysis for a broad range of
models. However, there are some limitations to be aware of:

* Adjoint sensitivity analysis is not fully supported for models with
  discontinuities (e.g. events).
  See `this issue <https://github.com/AMICI-dev/AMICI/issues/18>`__.

* Sensitivities may be incorrect for models with cascading events, where
  secondary events have a bolus.

* Sensitivities w.r.t. event trigger parameters may be incorrect for models
  with simultaneous events.

When working with a new model, it is recommended to check the sensitivities
against finite differences. This is demonstrated in the example notebooks.
This applies in particular to models with special features such as
events or to :term:`DAE` models.

Steady states
=============

AMICI provides various ways for computing steady states and steady state
sensitivities. However, there are some limitations to be aware of:

* A steady state is defined as a state where the time-derivative of the
  state variables is – and remains – zero. However, AMICI's steady state
  computation will stop as soon as the state time-derivative is smaller than
  the specified tolerance. This may lead to incorrect results
  for models with discontinuities (e.g. events) if the state time-derivative
  crosses the threshold before the discontinuity is reached.
  In this case, specify some large last output time instead of infinity to
  ensure that simulation does not stop prematurely.
  Note that AMICI only checks the model state (`x`), but not observables (`y`)
  or other model expressions (`w`) for steady state conditions.

* When computing steady states using Newton's method, only initial events
  will be handled. Subsequent events will be ignored.

Discontinuities
===============

AMICI is generally able to handle discontinuities in the model, such as events
or Boolean expressions. However, there are some limitations to be aware of:

* Not all event triggers are fully supported by the currently employed root
  finding algorithm. In particular, ``>`` and ``>=`` as well as ``<`` and
  ``<=`` cannot be distinguished. Also complex logical expressions where
  multiple operators are time-dependent may not be handled correctly.
  See, for example, `this issue <https://github.com/AMICI-dev/AMICI/issues/2707>`__.

* Discontinuities at the initial time point may not be fully supported.
  See, for example, `this issue <https://github.com/AMICI-dev/AMICI/issues/2724>`__.

Always check whether the simulation results are as expected.

Also, note the points above regarding sensitivity analysis and steady state
computation for models with discontinuities.
