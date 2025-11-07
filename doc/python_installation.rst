.. _amici_python_installation:

Installing the AMICI Python package
===================================

Short guide
+++++++++++

Installation of the AMICI Python package has the following prerequisites:

* Python>=3.11
* a C++20 compatible C++ compiler and a C compiler
  (e.g., g++>=10.1, clang>=13, Intel C++ compiler, mingw)

If these requirements are fulfilled and all relevant paths are setup properly,
AMICI can be installed using:

.. code-block:: bash

   pip3 install amici

If this worked, you can now import the Python module via::

   import amici

If this does not work for you, please follow the full instructions below.

.. note::

  To re-install a previously installed AMICI version with different
  build options or changed system libraries, pass the ``--no-cache-dir``
  option to ``pip`` to ensure a clean re-installation:

  .. code-block:: bash

     pip3 install --no-cache-dir amici


Installation on Linux
+++++++++++++++++++++

Ubuntu 22.04 / 24.04
--------------------

Install the AMICI dependencies via ``apt``
(this requires superuser privileges):

.. code-block:: bash

   sudo apt install python3-dev

   # optionally for HDF5 support:
   sudo apt install libhdf5-serial-dev

   # optionally for boost support (thread-specific CPU times, extended math functions, serialization)
   sudo apt install libboost-chrono-dev libboost-math-dev libboost-serialization-dev

Install AMICI:

.. code-block:: bash

   pip3 install amici


Arch Linux
----------

Install the AMICI dependencies via ``pacman``
(this requires superuser privileges):

.. code-block:: bash

   sudo pacman -S python gcc hdf5 boost-libs

Install AMICI:

.. code-block:: bash

   pip3 install amici

Alternatively:

1. Check if packages are already installed with the required versions for AMICI installation.

.. code-block:: bash

   sudo pacman -Si python gcc hdf5 boost-libs

2. Upgrade installed packages if required minimum versions are not satisfied for AMICI installation.

.. code-block:: bash

   sudo pacman -Su python gcc hdf5 boost-libs

4. Install AMICI:

.. code-block:: bash

   pip3 install amici


Installation on OSX
+++++++++++++++++++

Install the AMICI dependencies using homebrew:

.. code-block:: bash

    # optionally for HDF5 support:
    brew install hdf5

    # optionally for parallel simulations:
    brew install libomp
    # followed by either `brew link openmp` once,
    # or `export OpenMP_ROOT=$(brew --prefix)/opt/libomp"` where `OpenMP_ROOT` will have to be set during every re-installation of AMICI or any new model import

    # optionally for boost support (thread-specific CPU times, extended math functions, serialization)
    brew install boost && export BOOST_ROOT=$(brew --prefix)/opt/boost
    # followed by either `brew link boost` once,
    # or `export BOOST_ROOT=$(brew --prefix)/opt/boost"` where `BOOST_ROOT` will have to be set during every re-installation of AMICI or any new model import

Install AMICI:

.. code-block:: bash

    pip3 install amici


Installation on Windows
+++++++++++++++++++++++

Some general remarks:

* Consider using the `Windows Subsystem for Linux (WSL) <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`__ and follow the instructions for
  installation on linux.
* Install all libraries in a path not containing white spaces,
  e.g. directly under C:.
* Replace the following paths according to your installation.
* Slashes can be preferable to backslashes for some environment
  variables.
* See also [#425](https://github.com/AMICI-dev/amici/issues/425) for
  further discussion.

Using the Microsoft Visual Studio
---------------------------------

We assume that Visual Studio (not to be confused with Visual Studio Code)
is already installed. Using Visual Studio Installer, the following components
need to be included:

* Microsoft Visual C++ (MSVC).
  This is part of multiple packages, including Desktop Development with C++.
* Windows Universal C Runtime.
  This is an individual component and installs some DLLs that we need.

Further topics
++++++++++++++

Installation of development versions
------------------------------------

To install development versions which have not been released to PyPI yet,
you can install AMICI with ``pip`` directly from GitHub using:

.. code-block:: bash

    pip3 install -e git+https://github.com/AMICI-dev/amici.git@main#egg=amici\&subdirectory=python/sdist

Replace ``main`` by the branch or commit you want to install.

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
   dependencies listed in `setup.py <https://github.com/AMICI-dev/AMICI/blob/main/python/sdist/setup.py>`_
   manually, and try again. (This is because ``pip`` ``--install-option`` is
   applied to *all* installed packages, including dependencies.)


.. _amici_python_install_env_vars:

Custom installation
-------------------

Installation of the AMICI Python package can be customized using a number of
environment variables:

.. list-table:: Environment Variables for Custom Installation
   :header-rows: 1

   * - Variable
     - Purpose
     - Example
   * - ``SWIG``
     - Path to the :term:`SWIG` executable
     - ``SWIG=$HOME/bin/swig4.0``
   * - ``CC``
     - Setting the C(++) compiler
     - ``CC=/usr/bin/g++``
   * - ``CFLAGS``
     - Extra compiler flags used in every compiler invocation
     -
   * - ``AMICI_BLAS_USE_SCIPY_OPENBLAS``
     - Toggle using the OpenBLAS library provided by `scipy-openblas64`, on by default
       (if installed).
     - ``AMICI_BLAS_USE_SCIPY_OPENBLAS=FALSE``
   * - ``BLAS_CFLAGS``
     - Compiler flags for, e.g., BLAS include directories when using a non-default BLAS
     -
   * - ``BLAS_LIBS``
     - Flags for linking a non-default BLAS
     -
   * - ``ENABLE_GCOV_COVERAGE``
     - Set to build AMICI to generate code coverage information
     - ``ENABLE_GCOV_COVERAGE=TRUE``
   * - ``ENABLE_AMICI_DEBUGGING``
     - Set to build AMICI with debugging symbols
     - ``ENABLE_AMICI_DEBUGGING=TRUE``
   * - ``AMICI_PARALLEL_COMPILE``
     - Set to the number of parallel processes to be used for C(++) compilation (defaults to 1)
     - ``AMICI_PARALLEL_COMPILE=4``
   * - ``AMICI_TRY_ENABLE_HDF5``
     - Whether to build AMICI with HDF5-support if possible (Default: ``ON``)
     - ``AMICI_TRY_ENABLE_HDF5=OFF``


Installation under conda
------------------------

There is no amici conda recipe available yet. However, you can install AMICI
using pip in a conda environment.

Create a minimal conda environment via:

.. code-block:: bash

   conda create --name ENV_NAME pip python

Here, replace ``ENV_NAME`` by some name for the environment.

To activate the environment, run:

.. code-block:: bash

   source activate ENV_NAME

(and ``conda deactivate`` later to deactivate it again).

To install AMICI, now run:

.. code-block:: bash

   pip install amici

The ``pip`` option ``--no-cache`` may be helpful here to make sure the
installation is done completely anew.

Now, you are ready to use AMICI in the virtual environment.

Known issues:

* ``CMAKE_AR-NOTFOUND: not found``: Try ``conda install binutils``.

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
