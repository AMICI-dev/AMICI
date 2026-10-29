.. _versioning_policy:

Versioning policy
=================

Versioning
----------

We use `Semantic Versioning <http://semver.org/>`_ with the modifications
described under :ref:`deprecation_policy`.

.. _deprecation_policy:

Deprecation policy
------------------

AMICI aims to provide a stable API for users. However, not all features can be
maintained indefinitely. We will deprecate features in minor releases and
where possible, issue a warning when they are used. We will keep deprecated
features for at least six months after the release that includes the
respective deprecation warning and then remove them earliest in the next minor
or major release. If a deprecated feature is the source of a major bug, we may
remove it earlier.

Python compatibility
--------------------

We follow `SPEC 0 <https://scientific-python.org/specs/spec-0000/>`_
for the minimum supported versions of Python and of core dependencies such as
NumPy, SciPy, pandas, and Matplotlib.
We may keep supporting a Python version beyond the SPEC 0 schedule until
NumPy drops it.
