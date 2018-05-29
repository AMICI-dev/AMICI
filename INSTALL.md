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

### MATLAB

To use AMICI from MATLAB, start MATLAB and add the AMICI/matlab direcory to the MATLAB path. To add all toolbox directories to the MATLAB path, execute the matlab script

    installAMICI.m

To store the installation for further MATLAB session, the path can be saved via

    savepath

For the compilation of .mex files, MATLAB needs to be configured with a working C compiler. The C compiler needs to be installed and configured via:

    mex -setup c

For a list of supported compilers we refer to the mathworks documentation: [mathworks.com](http://mathworks.com/support/compilers/R2017a/index.html)

###  Python

To use AMICI from Python, install the module using pip  

    pip3 install amici
    
You can now import it as python module:

    import amici

### C++

To use AMICI from C++, run the  

    ./scripts/buildSundials.sh
    ./scripts/buildSuitesparse.sh
    ./scripts/buildAmici.sh
    
script to compile amici libary. The static library file can then be linked from 

    ./build/libamici.a
 
## Dependencies

The MATLAB interface requires the Mathworks Symbolic Toolbox for model generation via `amiwrap(...)`, but not for execution of precompiled models.

The Python and C++ interfaces require `cmake` and `cblas` to be installed.

The tools SUNDIALS and SuiteSparse shipped with AMICI do __not__ require explicit installation.

AMICI uses the following packages from SUNDIALS:

__CVODES__: the sensitivity-enabled ODE solver in SUNDIALS. Radu Serban and Alan C. Hindmarsh. _ASME 2005 International Design Engineering Technical Conferences and Computers and Information in Engineering Conference._ American Society of Mechanical Engineers, 2005. [PDF](http://proceedings.asmedigitalcollection.asme.org/proceeding.aspx?articleid=1588657)

__IDAS__

AMICI uses the following packages from SuiteSparse:

__Algorithm 907: KLU__, A Direct Sparse Solver for Circuit Simulation Problems. Timothy A. Davis, Ekanathan Palamadai Natarajan, _ACM Transactions on Mathematical Software_, Vol 37, Issue 6, 2010, pp 36:1 - 36:17. [PDF](http://dl.acm.org/authorize?305534)

__Algorithm 837: AMD__, an approximate minimum degree ordering algorithm, Patrick R. Amestoy, Timothy A. Davis, Iain S. Duff, _ACM Transactions on Mathematical Software_, Vol 30, Issue 3, 2004, pp 381 - 388. [PDF](http://dl.acm.org/authorize?733169)

__Algorithm 836: COLAMD__, a column approximate minimum degree ordering algorithm, Timothy A. Davis, John R. Gilbert, Stefan I. Larimore, Esmond G. Ng
_ACM Transactions on Mathematical Software_, Vol 30, Issue 3, 2004, pp 377 - 380. [PDF](http://dl.acm.org/authorize?734450)
