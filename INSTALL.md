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

To use AMICI from python, install the module and all other requirements using pip from pypi:

    pip3 install amici
    
You can now import it as python module:

    import amici

#### Anaconda

To use an Anaconda installation of python, proceeed as follows:

Since Anaconda provides own versions of packages which might not work with amici, create a minimum
virtual environment via:

    conda create --name ENV_NAME pip python

To activate the environment, do:

    source activate ENV_NAME

(and `conda deactivate` later to deactivate it again).

SWIG must be installed, as well as CBLAS. You can also use conda to install the latter locally, using:

    conda install -c conda-forge openblas

To install amici, now do:

    pip install amici

The option `--no-cache` may be helpful here to make sure the installation is done completely anew.

Now, you are ready to use amici in the virtual environment.

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
