.. _amici_python_faq:

FAQ
===

**Q**: I am trying to install the AMICI Python package, but installation
fails with something like

::

    error: no member named 'fill' in namespace 'std::ranges'

**A**: Try a more recent compiler. This error is caused by an older compiler
not supporting C++20 features used in the AMICI codebase.
