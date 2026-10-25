.. _amici_python_faq:

FAQ
===

**Q**: I am trying to install the AMICI Python package, but installation
fails with something like

::

    error: no member named 'fill' in namespace 'std::ranges'

**A**: Try a more recent compiler. This error is caused by an older compiler
not supporting C++20 features used in the AMICI codebase.

--------------

**Q**: Importing my model fails with something like
``ImportError: _someModelName.cpython-37m-x86_64-linux-gnu.so: undefined symbol: omp_get_thread_num``.

**A**: You probably installed the AMICI package with OpenMP support, but
did not have the relevant compiler/linker flags set when
importing/building the model. See :ref:`here <amici_python_openmp>`.
