# Installation

## Table of Contents
1. [Availability](#availability)
2. [Python](#python)
3. [MATLAB](#matlab)
4. [C++ only](#cpp)
5. [Dependencies](#dependencies)


<a name="availability"></a>
## Availability

The sources for AMICI are available as
- Source  [tarball](https://github.com/ICB-DCM/AMICI/tarball/master)
- Source  [zip](https://github.com/ICB-DCM/AMICI/zipball/master)
- GIT repository on  [github](https://github.com/ICB-DCM/AMICI)

A Python package is available on pypi, see below.

If AMICI was downloaded as a zip, it needs to be unpacked in a
convenient directory. If AMICI was obtained via cloning of the git
repository, no further unpacking is necessary.

### Obtaining AMICI via the GIT version control system
In order to always stay up-to-date with the latest AMICI versions,
simply pull it from our GIT repository and recompile it when a new
release is available. For more information about GIT checkout their
[website](http://git-scm.com/)

The GIT repository can currently be found at 
[https://github.com/ICB-DCM/AMICI](https://github.com/ICB-DCM/AMICI)
and a direct clone is possible via

    git clone https://github.com/ICB-DCM/AMICI.git AMICI


<a name="python"></a>
## Python

To use AMICI from python, install the module and all other requirements
using pip:

    pip3 install amici

You can now import it as python module:

    import amici

For cases where this installation fails, check below for special setups
and custom installations.
For Python-AMICI usage see 
[https://github.com/ICB-DCM/AMICI/blob/master/documentation/PYTHON.md](https://github.com/ICB-DCM/AMICI/blob/master/documentation/PYTHON.md).


### Installation of development versions

To install development versions which have not been released to pypi yet,
you can install AMICI with pip directly from GitHub using:

    pip3 install -e git+https://github.com/icb-dcm/amici.git@develop#egg=amici\&subdirectory=python/sdist

Replace `develop` by the branch or commit you want to install.

Note that this will probably not work on Windows which does not support
symlinks by default
(https://stackoverflow.com/questions/5917249/git-symlinks-in-windows/49913019#49913019).

### Light installation 

In case you only want to use the AMICI Python package for generating model code
for use with Matlab or C++ and don't want to bothered with any unnecessary 
dependencies, you can run

    pip3 install --install-option --no-clibs amici

Note, however, that you will not be able to compile any model into a 
Python extension with this installation.

NOTE: If you run into an error with above installation command, install
all AMICI dependencies listed in 
[`setup.py`](https://github.com/ICB-DCM/AMICI/blob/master/python/sdist/setup.py)
manually, and try again. (This is because `pip` `--install-option`s are
applied to *all* installed packages, including dependencies.)

### Anaconda

To use an Anaconda installation of python 
([https://www.anaconda.com/distribution/](https://www.anaconda.com/distribution/),
Python>=3.6), proceed as follows:

Since Anaconda provides own versions of some packages which might not
work with amici (in particular the gcc compiler), create a minimal
virtual environment via:

    conda create --name ENV_NAME pip python

Here, replace ENV_NAME by some name for the environment. To activate the
environment, do:

    source activate ENV_NAME

(and `conda deactivate` later to deactivate it again).

SWIG must be installed and available in your `PATH`, and a
CBLAS-compatible BLAS must be available. You can also use conda to
install the latter locally, using:

    conda install -c conda-forge openblas

To install AMICI, now do:

    pip install amici

The option `--no-cache` may be helpful here to make sure the
installation is done completely anew.

Now, you are ready to use AMICI in the virtual environment.

#### Anaconda on Mac

If the above installation does not work for you, try installing AMICI
via:

    CFLAGS="-stdlib=libc++" CC=clang CXX=clang pip3 install --verbose amici

This will use the `clang` compiler.

You will have to pass the same options when compiling any model later
on. This can be done by inserting the following code before calling
`sbml2amici`:
    
    import os
    os.environ['CC'] = 'clang'
    os.environ['CXX'] = 'clang'
    os.environ['CFLAGS'] = '-stdlib=libc++'

(For further discussion see https://github.com/ICB-DCM/AMICI/issues/357)


### Windows

To install AMICI on Windows using python, you can proceed as follows:

Some general remarks:

* Install all libraries in a path not containing white spaces,
  e.g. directly under C:.
* Replace the following paths according to your installation.
* Slashes can be preferable to backslashes for some environment
  variables.
* See also [#425](https://github.com/icb-dcm/amici/issues/425) for
  further discussion.

Then, follow these steps:

* A python environment for Windows is required. We recommend
  [Anaconda](https://www.anaconda.com/distribution/) with python >=3.6.
* Install [mingw64](https://sourceforge.net/projects/mingw-w64/files/latest/download)
  (32bit will succeed to compile, but fail during linking).
  During installation, select Version=8.1.0, Architecture=x64_64.
  Add the following directory to `PATH`:
    + `C:\mingw-w64\x86_64-8.1.0-posix-sjlj-rt_v6-rev0\mingw64\bin`
* Make sure that this is the compiler that is found by the system
  (e.g. `where gcc` in a `cmd` should point to this installation).
* Download CBLAS headers and libraries, e.g.
  [OpenBLAS](https://sourceforge.net/projects/openblas/files/v0.2.19/),
  binary distribution 0.2.19. Set the following environment variables:
    + `BLAS_CFLAGS=-IC:/OpenBLAS-v0.2.19-Win64-int32/include`
    + `BLAS_LIBS=-Wl,-Bstatic -LC:/OpenBLAS-v0.2.19-Win64-int32/lib -lopenblas -Wl,-Bdynamic`
* Install [SWIG](http://www.swig.org/download.html)
  (version swigwin-3.0.12 worked) and add the following directory to 
  `PATH`:
    + `C:\swigwin-3.0.12`
* Install AMICI using:

    `pip install --global-option="build_clib" --global-option="--compiler=mingw32" --global-option="build_ext" --global-option="--compiler=mingw32" amici --no-cache-dir --verbose`

Possible sources of errors:

* On recent Windows versions,
  `anaconda3\Lib\distutils\cygwinccompiler.py` fails linking
  `msvcr140.dll` with 
  `[...] x86_64-w64-mingw32/bin/ld.exe: cannot find -lmsvcr140`.
  This is not required for amici, so in `cygwinccompiler.py` 
  `return ['msvcr140']` can be changed to `return []`.
* If you use a python version where 
  [python/cpython#880](https://github.com/python/cpython/pull/880)
  has not been fixed yet, you need to disable 
  `define hypot _hypot in anaconda3\include/pyconfig.h` yourself.
* `import amici` in python resulting in the very informative 
  
  > ImportError: DLL load failed: The specified module could not be found.
  
  means that some amici module dependencies were not found (not the
  AMICI module itself).
  [DependencyWalker](http://www.dependencywalker.com/) will show you
  which ones.
  
  Support for msvc is experimental. [installOpenBLAS.ps1](https://github.com/ICB-DCM/AMICI/blob/master/scripts/installOpenBLAS.ps1) and [compileBLAS.cmd](https://github.com/ICB-DCM/AMICI/blob/master/scripts/compileBLAS.cmd) may serve as guidance on how to install openBLAS using msvc.

### Custom installation

AMICI Python package installation can be customized using a number of
environment variables:

|Variable | Purpose | Example |
|---|---|---|
|`CC`| Setting the C(++) compiler | `CC=/usr/bin/g++`| 
|`CFLAGS`| Extra compiler flags used in every compiler call | | 
|`BLAS_CFLAGS`| Compiler flags for, e.g. BLAS include directories | | 
|`BLAS_LIBS`| Flags for linking BLAS | | 
|`ENABLE_GCOV_COVERAGE`| Set to build AMICI to provide code coverage information | `ENABLE_GCOV_COVERAGE=TRUE`| 
|`ENABLE_AMICI_DEBUGGING`| Set to build AMICI with debugging symbols | `ENABLE_AMICI_DEBUGGING=TRUE`| 
|`AMICI_PARALLEL_COMPILE`| Set to the number of parallel processes to be used for C(++) file compilation (defaults to 1)| `AMICI_PARALLEL_COMPILE=4`|


<a name="matlab"></a>
## MATLAB

To use AMICI from MATLAB, start MATLAB and add the `AMICI/matlab`
directory to the MATLAB path. To add all toolbox directories to the
MATLAB path, execute the matlab script

    installAMICI.m

To store the installation for further MATLAB session, the path can be
saved via

    savepath

For the compilation of .mex files, MATLAB needs to be configured with a
working C++ compiler. The C++ compiler needs to be installed and
configured via:

    mex -setup c++

For a list of supported compilers we refer to the mathworks 
documentation: 
[mathworks.com](http://mathworks.com/support/compilers/R2018b/index.html)
Note that Microsoft Visual Studio compilers are currently not supported.


<a name="cpp"></a>
## C++ only

To use AMICI from C++, run the

    ./scripts/buildSundials.sh
    ./scripts/buildSuitesparse.sh
    ./scripts/buildAmici.sh

script to compile AMICI library.

**NOTE**: On some systems, the CMake executable may be named something
other than `cmake`. In this case, set the `CMAKE` environment variable
to the correct name (e.g. `export CMAKE=cmake3`, in case you have CMake
available as `cmake3`).

The static library file can then be linked from

    ./build/libamici.a

In CMake-based packages, amici can be linked via

    find_package(Amici)

### Optional SuperLU_MT support

To build AMICI with SuperLU_MT support, run

    ./scripts/buildSuperLUMT.sh
    ./scripts/buildSundials.sh
    cd build/
    cmake -DSUNDIALS_SUPERLUMT_ENABLE=ON ..
    make


<a name="dependencies"></a>
## Dependencies

### General

The tools SUNDIALS and SuiteSparse shipped with AMICI do __not__ require
explicit installation.

AMICI uses the following packages from SUNDIALS:

__CVODES__: the sensitivity-enabled ODE solver in SUNDIALS. Radu Serban
and Alan C. Hindmarsh. _ASME 2005 International Design Engineering
Technical Conferences and Computers and Information in Engineering
Conference._ American Society of Mechanical Engineers, 2005. 
[PDF](http://proceedings.asmedigitalcollection.asme.org/proceeding.aspx?articleid=1588657)

__IDAS__

AMICI uses the following packages from SuiteSparse:

__Algorithm 907: KLU__, A Direct Sparse Solver for Circuit Simulation
Problems. Timothy A. Davis, Ekanathan Palamadai Natarajan, 
_ACM Transactions on Mathematical Software_, Vol 37, Issue 6, 2010,
pp 36:1 - 36:17. [PDF](http://dl.acm.org/authorize?305534)

__Algorithm 837: AMD__, an approximate minimum degree ordering
algorithm, Patrick R. Amestoy, Timothy A. Davis, Iain S. Duff,
_ACM Transactions on Mathematical Software_, Vol 30, Issue 3, 2004,
pp 381 - 388. [PDF](http://dl.acm.org/authorize?733169)

__Algorithm 836: COLAMD__, a column approximate minimum degree ordering
algorithm, Timothy A. Davis, John R. Gilbert, Stefan I. Larimore,
Esmond G. Ng _ACM Transactions on Mathematical Software_, Vol 30,
Issue 3, 2004, pp 377 - 380. [PDF](http://dl.acm.org/authorize?734450)

#### libsbml

To import Systems Biology Markup Language ([SBML](http://sbml.org/))
models, AMICI relies on the Python or MATLAB SBML library.

#### Math Kernel Library (MKL)

The python and C++ interfaces require a system installation of a `BLAS`.
AMICI has been tested with various native and general purpose MKL
implementations such as Accelerate, Intel MKL, cblas, openblas, atlas.
The matlab interface uses the MATLAB MKL, which requires no separate
installation.

On Ubuntu, this requirement can be satisfied with

    apt install libatlas-base-dev

On Fedora (32):

    sudo dnf install blas-devel

#### C++ compiler

All AMICI installations require a C++11-compatible C++ compiler.
AMICI has been tested with g++, mingw, clang and the Intel compiler.
Visual C++ is not officially supported, but may work.

#### HDF5

The python and C++ interfaces provide routines to read and write options
and results in [hdf5](https://support.hdfgroup.org/HDF5/) format. 
For the python interface, the installation of hdf5 is optional, but for
the C++ interace it is currently required.

HDF5 can be installed using package managers such as
[brew](https://brew.sh) or [apt](https://wiki.debian.org/Apt):

    brew install hdf5

or

    apt-get install libhdf5-serial-dev

#### SWIG

The python interface requires [SWIG](http://www.swig.org), which has to
be installed by the user. As root user, SWIG can be installed using
package managers such as [brew](https://brew.sh) or
[apt](https://wiki.debian.org/Apt):

    brew install swig

or

    apt-get install swig3.0

Or by non-root users, using `scripts/downloadAndBuildSwig.sh` from the
AMICI repository (not included in the PyPI package). The binary
directory has to be added to the `PATH` environment variable, or `SWIG`
has to be set as described in the following section.


##### Using a non-default SWIG executable

We note here that some linux package managers may provide swig
executables as `swig3.0`, but installation as `swig` is required. This
can be fixed as root user using, e.g., symbolic links:

    mkdir -p ~/bin/ && ln -s $(which swig3.0) ~/bin/swig && export PATH=~/bin/:$PATH

Non-root users can set the `SWIG` environment variable to the full
path of the desired SWIG executable. This variable has be set during
AMICI package installation as well as during model compilation.

### Matlab

The MATLAB interface requires the Mathworks Symbolic Toolbox for model
generation via `amiwrap(...)`, but not for execution of precompiled
models. Currently MATLAB R2018a or newer is not supported (see 
[https://github.com/ICB-DCM/AMICI/issues/307](https://github.com/ICB-DCM/AMICI/issues/307)).

The Symbolic Toolbox requirement can be circumvented by performing model
import using the Python interface. The result code can then be used from
Matlab. 


### Python

The python interface requires python 3.6 or newer and a cblas-compatible
BLAS library to be installed. Windows installations via pip are
currently not supported, but users may try to install amici using the
build scripts provided for the C++ interface (these will by default
automatically install the python module).

The python interface depends on some additional packages, e.g. `numpy`.
They are automatically installed when installing the python package.

### C++

The C++ interface requires `cmake` and a cblas-compatible BLAS to be
installed.


### Optional

__SuperLU_MT__, "a general purpose library for the direct solution of large,
sparse, nonsymmetric systems of linear equations"
(https://crd-legacy.lbl.gov/~xiaoye/SuperLU/#superlu_mt).
SuperLU_MT is optional and is so far only available from the C++ interface.
