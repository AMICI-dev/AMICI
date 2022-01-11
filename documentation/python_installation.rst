.. _amici_python_installation:

Installing the AMICI Python package
===================================

Short guide
+++++++++++

Installation of the AMICI Python package has the following prerequisites:

* Python>=3.8
* :term:`SWIG`>=3.0
* CBLAS compatible BLAS library
  (e.g., OpenBLAS, CBLAS, Atlas, Accelerate, Intel MKL)
* a C++14 compatible C++ compiler and a C compiler
  (e.g., g++, clang, Intel C++ compiler, mingw)

If these requirements are fulfilled and all relevant paths are setup properly,
AMICI can be installed using:

.. code-block:: bash

   pip3 install amici

If this worked, you can now import the Python module via::

   import amici

If this does not work for you, please follow the full instructions below.

Installation on Linux
+++++++++++++++++++++

Ubuntu 20.04
------------

Install the AMICI dependencies via ``apt``
(this requires superuser privileges):

.. code-block:: bash

   sudo apt install libatlas-base-dev swig

   # optionally for HDF5 support:
   sudo apt install libhdf5-serial-dev

Install AMICI:

.. code-block:: bash

   pip3 install amici

Fedora 32
---------

Install the AMICI dependencies via ``apt``
(this requires superuser privileges):

.. code-block:: bash

   sudo dnf install blas-devel swig

Install AMICI:

.. code-block:: bash

   pip3 install amici

Installation on OSX
+++++++++++++++++++

Install the AMICI dependencies using homebrew:

.. code-block:: bash

    brew install swig

    # optionally for HDF5 support:
    brew install hdf5

    # optionally for parallel simulations:
    brew install libomp

Install AMICI:

.. code-block:: bash

    pip3 install amici


Installation on Windows
+++++++++++++++++++++++

Some general remarks:

* Install all libraries in a path not containing white spaces,
  e.g. directly under C:.
* Replace the following paths according to your installation.
* Slashes can be preferable to backslashes for some environment
  variables.
* See also [#425](https://github.com/AMICI-dev/amici/issues/425) for
  further discussion.

Using the MinGW compilers
-------------------------

* Install `MinGW-W64 <https://sourceforge.net/projects/mingw-w64/files/>`_
  (the 32bit version will succeed to compile, but fail during linking).

  MinGW-W64 GCC-8.1.0 for ``x86_64-posix-sjlj``
  (`direct link <https://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win64/Personal%20Builds/mingw-builds/8.1.0/threads-posix/sjlj/x86_64-8.1.0-release-posix-sjlj-rt_v6-rev0.7z/download>`_) has been shown to work on Windows 7 and 10 test systems.

* Add the following directory to your ``PATH``:
  ``C:\mingw-w64\x86_64-8.1.0-posix-sjlj-rt_v6-rev0\mingw64\bin``

* Make sure that this is the compiler that is found by the system
  (e.g. ``where gcc`` in a ``cmd`` should point to this installation).

* Download CBLAS headers and libraries, e.g.
  `OpenBLAS <https://sourceforge.net/projects/openblas/files/v0.2.19/>`_,
  binary distribution 0.2.19.

  Set the following environment variables:

  + ``BLAS_CFLAGS=-IC:/OpenBLAS-v0.2.19-Win64-int32/include``
  + ``BLAS_LIBS=-Wl,-Bstatic -LC:/OpenBLAS-v0.2.19-Win64-int32/lib -lopenblas -Wl,-Bdynamic``

* Install `SWIG <http://www.swig.org/download.html>`_
  and add the SWIG directory to ``PATH``
  (e.g. ``C:\swigwin-3.0.12`` for version 3.0.12)

* Install AMICI using:

  .. code-block:: bash

     pip install --global-option="build_clib" \
                 --global-option="--compiler=mingw32" \
                 --global-option="build_ext" \
                 --global-option="--compiler=mingw32" \
                 amici --no-cache-dir --verbose`

.. note::

   **Possible sources of errors:**

   * On recent Windows versions,
     ``anaconda3\Lib\distutils\cygwinccompiler.py`` fails linking
     ``msvcr140.dll`` with
     ``[...] x86_64-w64-mingw32/bin/ld.exe: cannot find -lmsvcr140``.
     This is not required for amici, so in ``cygwinccompiler.py``
     ``return ['msvcr140']`` can be changed to ``return []``.

   * If you use a python version where
     `python/cpython#880 <https://github.com/python/cpython/pull/880>`_
     has not been fixed yet, you need to disable
     ``define hypot _hypot`` in ``anaconda3\include/pyconfig.h`` yourself.

   * ``import amici`` in Python resulting in the very informative

       ImportError: DLL load failed: The specified module could not be found.

     means that some amici module dependencies were not found (not the
     AMICI module itself).
     `DependencyWalker <http://www.dependencywalker.com/>`_ can show you
     which ones.

Using the Microsoft Visual Studio
---------------------------------

.. note:: Support for MSVC is experimental.

We assume that Visual Studio (not to be confused with Visual Studio Code)
is already installed. Using Visual Studio Installer, the following components
need to be included:

* Microsoft Visual C++ (MSVC).
  This is part of multiple packages, including Desktop Development with C++.
* Windows Universal C Runtime.
  This is an individual component and installs some DLLs that we need.

OpenBLAS
^^^^^^^^

There are prebuilt OpenBLAS binaries available, but they did not seem to work
well here. Therefore, we recommend building OpenBLAS from scratch.

To build OpenBLAS, download the following scripts from the AMICI repository:

* https://github.com/AMICI-dev/AMICI/blob/master/scripts/installOpenBLAS.ps1
* https://github.com/AMICI-dev/AMICI/blob/master/scripts/compileBLAS.cmd

The first script needs to be called in Powershell, and it needs to call
``compileBLAS.cmd``, so you will need to modify line 11:

    C: \\Users\\travis\\build\\AMICI\\scripts\\compileBLAS.cmd

so that it matches your directory structure.
This will download OpenBLAS and compile it, creating

    C:\\BLAS\lib\\openblas.lib
    C:\\BLAS\\bin\\openblas.dll

You will also need to define two environment variables:

.. code-block:: text

   BLAS_LIBS="/LIBPATH:C:\BLAS\lib openblas.lib"
   BLAS_CFLAGS="/IC:/BLAS/OpenBLAS-0.3.12/OpenBLAS-0.3.12"

One way to do that is to run a PowerShell script with the following commands:

.. code-block:: text

   [System.Environment]::SetEnvironmentVariable("BLAS_LIBS", "/LIBPATH:C:/BLAS/lib openblas.lib", [System.EnvironmentVariableTarget]::User)
   [System.Environment]::SetEnvironmentVariable("BLAS_LIBS", "/LIBPATH:C:/BLAS/lib openblas.lib", [System.EnvironmentVariableTarget]::Process)
   [System.Environment]::SetEnvironmentVariable("BLAS_CFLAGS", "-IC:/BLAS/OpenBLAS-0.3.12/OpenBLAS-0.3.12", [System.EnvironmentVariableTarget]::User)
   [System.Environment]::SetEnvironmentVariable("BLAS_CFLAGS", "-IC:/BLAS/OpenBLAS-0.3.12/OpenBLAS-0.3.12", [System.EnvironmentVariableTarget]::Process)

The call ending in ``Process`` sets the environment variable in the current
process, and it is no longer in effect in the next process. The call ending in
``User`` is permanent, and takes effect the next time the user logs on.

Now you need to make sure that all required DLLs are within the scope of the
``PATH`` variable. In particular, the following directories need to be included
in ``PATH``:

    C:\\BLAS\\bin
    C:\\Program Files (x86)\\Windows Kits\\10\\Redist\\ucrt\\DLLs\\x64

The first one is needed for ``openblas.dll`` and the second is needed for the
Windows Universal C Runtime.

If any DLLs are missing in the ``PATH`` variable, Python will return the
following error upon ``import amici``:

    ImportError: DLL load failed: The specified module could not be found.

This can be tested using the "where" command. For example

    where openblas.dll

should return

    C:\\BLAS\\bin\\openblas.dll

Almost all of the DLLs are standard Windows DLLs and should be included in
either Windows or Visual Studio. But, in case it is necessary to test this,
here is a list of some DLLs required by AMICI (when compiled with MSVC):

* ``openblas.dll``
* ``python37.dll``
* ``MSVCP140.dll``
* ``KERNEL32.dll``
* ``VCRUNTIME140_1.dll``
* ``VCRUNTIME140.dll``
* ``api-ms-win-crt-convert-l1-1-0.dll``
* ``api-ms-win-crt-heap-l1-1-0.dll``
* ``api-ms-win-crt-stdio-l1-1-0.dll``
* ``api-ms-win-crt-string-l1-1-0.dll``
* ``api-ms-win-crt-runtime-l1-1-0.dll``
* ``api-ms-win-crt-time-l1-1-0.dll``
* ``api-ms-win-crt-math-l1-1-0.dll``

``MSVCP140.dll``, ``VCRUNTIME140.dll``, and ``VCRUNTIME140_1.dll`` are needed
by MSVC (see Visual Studio above). ``KERNEL32.dll`` is part of Windows and in
``C:\Windows\System32``. The ``api-ms-win-crt-XXX-l1-1-0.dll`` are needed by
``openblas.dll`` and are part of the Windows Universal C Runtime.

.. note::

    Since Python 3.8, the library directory needs to be set either from Python:

    .. code-block:: python

        import os
        # directory containing `openblas.dll`
        os.add_dll_directory("C:\\BLAS\\bin")
        import amici

    or via the environment variable ``AMICI_DLL_DIRS``.

Further topics
++++++++++++++

Installation of development versions
------------------------------------

To install development versions which have not been released to PyPI yet,
you can install AMICI with ``pip`` directly from GitHub using:

.. code-block:: bash

    pip3 install -e git+https://github.com/AMICI-dev/amici.git@develop#egg=amici\&subdirectory=python/sdist

Replace ``develop`` by the branch or commit you want to install.

Note that this will only work on Windows if you have enabled developer mode,
because symlinks are not supported by default
(`more information <https://stackoverflow.com/questions/5917249/git-symlinks-in-windows/49913019#49913019>`_).

Light installation
------------------

In case you only want to use the AMICI Python package for generating model code
for use with Matlab or C++ and don't want to bothered with any unnecessary
dependencies, you can run

.. code-block:: bash

   pip3 install --install-option --no-clibs amici

.. note::

   Following this installation, you will not be able to simulate the imported
   models in Python.

.. note::

   If you run into an error with above installation command, install all AMICI
   dependencies listed in `setup.py <https://github.com/AMICI-dev/AMICI/blob/master/python/sdist/setup.py>`_
   manually, and try again. (This is because ``pip`` ``--install-option`` is
   applied to *all* installed packages, including dependencies.)


Custom installation
-------------------

Installation of the AMICI Python package can be customized using a number of
environment variables:

+----------------------------+----------------------------------+---------------------------------+
| Variable                   | Purpose                          | Example                         |
+============================+==================================+=================================+
| ``SWIG``                   | Path to the :term:`SWIG`         | ``SWIG=$HOME/bin/swig4.0``      |
|                            | executable                       |                                 |
+----------------------------+----------------------------------+---------------------------------+
| ``CC``                     | Setting the C(++) compiler       | ``CC=/usr/bin/g++``             |
+----------------------------+----------------------------------+---------------------------------+
| ``CFLAGS``                 | Extra compiler flags used in     |                                 |
|                            | every compiler invocation        |                                 |
+----------------------------+----------------------------------+---------------------------------+
| ``BLAS_CFLAGS``            | Compiler flags for, e.g. BLAS    |                                 |
|                            |  include directories             |                                 |
+----------------------------+----------------------------------+---------------------------------+
| ``BLAS_LIBS``              | Flags for linking BLAS           |                                 |
+----------------------------+----------------------------------+---------------------------------+
| ``ENABLE_GCOV_COVERAGE``   | Set to build AMICI to generate   | ``ENABLE_GCOV_COVERAGE=TRUE``   |
|                            | code coverage information        |                                 |
+----------------------------+----------------------------------+---------------------------------+
| ``ENABLE_AMICI_DEBUGGING`` | Set to build AMICI with          | ``ENABLE_AMICI_DEBUGGING=TRUE`` |
|                            | debugging symbols                |                                 |
+----------------------------+----------------------------------+---------------------------------+
| ``AMICI_PARALLEL_COMPILE`` | Set to the number of parallel    | ``AMICI_PARALLEL_COMPILE=4``    |
|                            | processes to be used for C(++)   |                                 |
|                            | compilation (defaults to 1)      |                                 |
+----------------------------+----------------------------------+---------------------------------+

Installation under Anaconda
---------------------------

To use an Anaconda installation of Python
`https://www.anaconda.com/distribution/ <https://www.anaconda.com/distribution/>`_,
Python>=3.7), proceed as follows:

Since Anaconda provides own versions of some packages which might not
work with AMICI (in particular the ``gcc`` compiler), create a minimal
virtual environment via:

.. code-block:: bash

   conda create --name ENV_NAME pip python

Here, replace ``ENV_NAME`` by some name for the environment.

To activate the environment, run:

.. code-block:: bash

   source activate ENV_NAME

(and ``conda deactivate`` later to deactivate it again).

:term:`SWIG` must be installed and available in your ``PATH``, and a
CBLAS-compatible BLAS must be available. You can also use conda to
install the latter locally, using:

.. code-block:: bash

   conda install -c conda-forge openblas

To make AMICI use openblas, set the following environment variable:

.. code-block:: bash

   export BLAS_LIBS=-lopenblas

``BLAS_LIBS`` needs to be set during installation of the AMICI package, as
well as during any future model import.

To install AMICI, now run:

.. code-block:: bash

   pip install amici

The ``pip`` option ``--no-cache`` may be helpful here to make sure the
installation is done completely anew.

Now, you are ready to use AMICI in the virtual environment.

.. note::

   **Anaconda on Mac**

   If the above installation does not work for you, try installing AMICI via:

   .. code-block:: bash

      CFLAGS="-stdlib=libc++" CC=clang CXX=clang pip3 install --verbose amici

   This will use the ``clang`` compiler.

   You will have to pass the same options when compiling any model later
   on. This can be done by inserting the following code before model import:

   .. code-block:: python

      import os
      os.environ['CC'] = 'clang'
      os.environ['CXX'] = 'clang'
      os.environ['CFLAGS'] = '-stdlib=libc++'

   (For further discussion see https://github.com/AMICI-dev/AMICI/issues/357)


Optional Boost support
----------------------

`Boost <https://www.boost.org/>`_ is an optional C++ dependency only required
for special functions (including e.g. gamma derivatives) in the Python
interface. Boost can be installed via package managers via

.. code-block:: bash

    apt-get install libboost-math-dev

or

.. code-block:: bash

    brew install boost

As only headers are required, also a
`source code <https://www.boost.org/doc/libs/1_66_0/more/getting_started/unix-variants.html>`_
download suffices. The compiler must be able to find the module in the search
path.
