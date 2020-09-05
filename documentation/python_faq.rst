FAQ
===

**Q**: I am trying to install the AMICI Python package, but installation
fails with something like

::

   amici/src/cblas.cpp:16:13: fatal error: cblas.h: No such file or directory
   #include <cblas.h>
            ^~~~~~~~~
   compilation terminated.
   error: command 'x86_64-linux-gnu-gcc' failed with exit status 1

**A**: You will have to install a CBLAS-compatible BLAS library and/or
set ``BLAS_CFLAGS`` as described in the installation guides.

--------------

**Q**: Importing my model fails with something like
``ImportError: _someModelName.cpython-37m-x86_64-linux-gnu.so: undefined symbol: omp_get_thread_num``.

**A**: You probably installed the AMICI package with OpenMP support, but
did not have the relevant compiler/linker flags set when
importing/building the model. See `Python-AMICI
guide <python_interface.rst#model-compilation>`__.
