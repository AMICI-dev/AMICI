# Installation

## Availability

The sources for AMICI are accessible as
- Source  [tarball](https://github.com/ICB-DCM/AMICI/tarball/master)
- Source  [zip](https://github.com/ICB-DCM/AMICI/zipball/master)
- GIT repository on  [github](https://github.com/ICB-DCM/AMICI)

### Obtaining AMICI via the GIT versioning system
In order to always stay up-to-date with the latest AMICI versions, simply pull it from our GIT repository and
recompile it when a new release is available. For more information about GIT checkout their [website](http://git-scm.com/)

The GIT repository can currently be found at https://github.com/ICB-DCM/AMICI and a direct clone is possible via

    git clone https://github.com/ICB-DCM/AMICI.git AMICI

## Installation

If AMICI was downloaded as a zip, it needs to be unpacked in a convenient directory. If AMICI was obtained via cloning of the git repository, no further unpacking is necessary.

### Dependencies

The MATLAB interface only depends on the symbolic toolbox, which is needed for model compilation, but not simulation.

#### Symbolic Engine

The MATLAB interface requires the symbolic toolbox for model compilation. The symbolic toolbox is not required for model simulation.

#### Math Kernel Library (MKL)

The python and C++ interfaces require a system installation of  `BLAS`. AMICI has been tested with various native and general purpose MKL implementations such as
Accelerate, Intel MKL, cblas, openblas, atlas. The matlab interface uses the MATLAB MKL, which requires no prior installation.

#### HDF5

The python and C++ interfaces provide routines to read and write options and results in [hdf5](https://support.hdfgroup.org/HDF5/) format. For the python interface, the installation of hdf5 is optional, but for the C++ interace it is required. 
HDF can be installed using package managers such as [brew](https://brew.sh) or [apt](https://wiki.debian.org/Apt):

    brew install hdf5

or

    apt-get install libhdf5-serial-dev


#### SWIG

The python interface requires [SWIG](http://www.swig.org), which has to be installed by the user. 
Swig can be installed using package managers such as [brew](https://brew.sh) or [apt](https://wiki.debian.org/Apt):

    brew install swig

or

    apt-get install swig3.0

We note here that some linux package managers may provide swig executables as `swig3.0`, but installation as `swig` is required. This can be fixed using, e.g., symbolic links:

    mkdir -p ~/bin/ && ln -s $(which swig3.0) ~/bin/swig && export PATH=~/bin/:$PATH

#### python packages

The python interface requires the python packages `pkgconfig` and `numpy` to be installed before AMICI can be installed. These can be installed via `pip`:

    pip3 install pkgconfig numpy

### MATLAB

To use AMICI from MATLAB, start MATLAB and add the AMICI/matlab direcory to the MATLAB path. To add all toolbox directories to the MATLAB path, execute the matlab script

    installAMICI.m

To store the installation for further MATLAB session, the path can be saved via

    savepath

For the compilation of .mex files, MATLAB needs to be configured with a working C compiler. The C compiler needs to be installed and configured via:

    mex -setup c

For a list of supported compilers we refer to the mathworks documentation: [mathworks.com](http://mathworks.com/support/compilers/R2018b/index.html)
Note that Microsoft Visual Studio compilers are currently not supported.

###  python

To use AMICI from python, install the module and all other requirements using pip:

    pip3 install amici

You can now import it as python module:

    import amici

In case you only want to use the AMICI Python package for generating model code
for use with Matlab or C++ and don't want to bothered with any unnecessary 
dependencies, you can run

    pip3 install --install-option --no-clibs amici

Note, however, that you will not be able to compile any model into a Python
extension with this installation.

#### Anaconda

To use an Anaconda installation of python (https://www.anaconda.com/distribution/, Python>=3.6), proceed as follows:

Since Anaconda provides own versions of some packages which might not work with amici (in particular the gcc compiler), create a minimal virtual environment via:

    conda create --name ENV_NAME pip python

Here, replace ENV_NAME by some name for the environment. To activate the environment, do:

    source activate ENV_NAME

(and `conda deactivate` later to deactivate it again).

SWIG must be installed and available in your `PATH`, and a CBLAS-compatible BLAS must be available. You can also use conda to install the latter locally, using:

    conda install -c conda-forge openblas

To install AMICI, now do:

    pip install amici

The option `--no-cache` may be helpful here to make sure the installation is done completely anew.

Now, you are ready to use AMICI in the virtual environment.

#### Windows

To install AMICI on Windows using python, you can proceed as follows:

Some general remarks:

* Install all libraries in a path not containing white spaces, e.g. directly under C:.
* Replace the following paths according to your installation.
* Slashes can be preferable to backslashes for some environment variables.
* See also [#425](https://github.com/icb-dcm/amici/issues/425) for further discussion.

Then, follow these steps:

* A python environment for Windows is required. We recommend [Anaconda](https://www.anaconda.com/distribution/) with python >=3.6.
* Install [mingw64](https://sourceforge.net/projects/mingw-w64/files/latest/download) (32bit will succeed to compile, but fail during linking). During installation, select Version=8.1.0, Architecture=x64_64. Add the following directory to `PATH`:
    + `C:\mingw-w64\x86_64-8.1.0-posix-sjlj-rt_v6-rev0\mingw64\bin`
* Make sure that this is the compiler that is found by the system (e.g. `where gcc` in a `cmd` should point to this installation).
* Download CBLAS headers and libraries, e.g. [OpenBLAS](https://sourceforge.net/projects/openblas/files/v0.2.19/), binary distribution 0.2.19. Set the following environment variables:
    + `BLAS_CFLAGS=-IC:/OpenBLAS-v0.2.19-Win64-int32/include`
    + `BLAS_LIBS=-Wl,-Bstatic -LC:/OpenBLAS-v0.2.19-Win64-int32/lib -lopenblas -Wl,-Bdynamic`
* Install [SWIG](http://www.swig.org/download.html) (version swigwin-3.0.12 worked) and add the following directory to `PATH`:
    + `C:\swigwin-3.0.12`
* Install AMICI using:

    `pip install --global-option="build_clib" --global-option="--compiler=mingw32" --global-option="build_ext" --global-option="--compiler=mingw32" amici --no-cache-dir --verbose`

Possible sources of errors:

* On recent Windows versions, `anaconda3\Lib\distutils\cygwinccompiler.py` fails linking `msvcr140.dll` with `[...] x86_64-w64-mingw32/bin/ld.exe: cannot find -lmsvcr140`. This is not required for amici, so in `cygwinccompiler.py` `return ['msvcr140']` can be changed to `return []`.
* If you use a python version where python/cpython#880 has not been fixed yet, you need to disable `#define hypot _hypot in anaconda3\include/pyconfig.h` yourself.
* `import amici` in python resulting in the very informative "ImportError: DLL load failed: The specified module could not be found." means that some amici module dependencies were not found (not the amici module itself). [DependencyWalker](http://www.dependencywalker.com/) will show you which ones.

### C++

To use AMICI from C++, run the

    ./scripts/buildSundials.sh
    ./scripts/buildSuitesparse.sh
    ./scripts/buildAmici.sh

script to compile amici libary. The static library file can then be linked from 

    ./build/libamici.a

In CMake-based packages, amici can be linked via

    find_package(Amici)

## Dependencies

The MATLAB interface requires the Mathworks Symbolic Toolbox for model generation via `amiwrap(...)`, but not for execution of precompiled models. Currently MATLAB R2018a or newer is not supported (see https://github.com/ICB-DCM/AMICI/issues/307)

The python interface requires python 3.6 or newer and `cblas` library to be installed. Windows installations via pip are currently not supported, but users may try to install amici using the build scripts provided for the C++ interface (these will by default automatically install the python module).

The C++ interface requires `cmake` and `cblas` to be installed.

The tools SUNDIALS and SuiteSparse shipped with AMICI do __not__ require explicit installation.

AMICI uses the following packages from SUNDIALS:

__CVODES__: the sensitivity-enabled ODE solver in SUNDIALS. Radu Serban and Alan C. Hindmarsh. _ASME 2005 International Design Engineering Technical Conferences and Computers and Information in Engineering Conference._ American Society of Mechanical Engineers, 2005. [PDF](http://proceedings.asmedigitalcollection.asme.org/proceeding.aspx?articleid=1588657)

__IDAS__

AMICI uses the following packages from SuiteSparse:

__Algorithm 907: KLU__, A Direct Sparse Solver for Circuit Simulation Problems. Timothy A. Davis, Ekanathan Palamadai Natarajan, _ACM Transactions on Mathematical Software_, Vol 37, Issue 6, 2010, pp 36:1 - 36:17. [PDF](http://dl.acm.org/authorize?305534)

__Algorithm 837: AMD__, an approximate minimum degree ordering algorithm, Patrick R. Amestoy, Timothy A. Davis, Iain S. Duff, _ACM Transactions on Mathematical Software_, Vol 30, Issue 3, 2004, pp 381 - 388. [PDF](http://dl.acm.org/authorize?733169)

__Algorithm 836: COLAMD__, a column approximate minimum degree ordering algorithm, Timothy A. Davis, John R. Gilbert, Stefan I. Larimore, Esmond G. Ng
_ACM Transactions on Mathematical Software_, Vol 30, Issue 3, 2004, pp 377 - 380. [PDF](http://dl.acm.org/authorize?734450)
