# C++ Interface

The [Python interface](https://amici.readthedocs.io/en/latest/PYTHON.html)
and @ref matlab_interface can translate the model
definition into C++ code, which is then compiled into a `.mex` file or a Python
module. Advanced users can also use this code within stand-alone C/C++
application for use in other environments (e.g. on high performance computing
systems). This section will give a short overview over the generated files and
provide a brief introduction of how this code can be included in other
applications. For more details see the function level documentation, e.g. in
http://icb-dcm.github.io/AMICI/amici_8cpp.html. 


## Generation of C++ model files

`amiwrap.m` (Matlab) and `amici.SbmlImporter.sbml2amici` (Python) import a
model and create C++ code for various model functions. The model
source files are written to `${AMICI_ROOT_DIR}/models/${MODEL_NAME}` by
default.

The content of a model source directory might look something like this (given
`MODEL_NAME=model_steadystate`): 

```
CMakeLists.txt
hashes.mat 
main.cpp 
model_steadystate_deltaqB.cpp 
model_steadystate_deltaqB.h 
[... many more files model_steadystate_*.(cpp|h|md5|o) ]
wrapfunctions.cpp 
wrapfunctions.h 
model_steadystate.h 
```

These files provide the implementation of a model-specific subclass of
`amici::Model`. The `CMakeLists.txt` file can be used to build the library.
`main.cpp` contains a simple scaffold for running a model simulation from C++.
See next section for more details on these files.


## Compiling and linking

In your application, you need to compile and link every model 
`${AMICI_ROOT_DIR}/models/${MODEL_NAME}/*.cpp`,
`${AMICI_ROOT_DIR}/src/*.cpp`, the SUNDIALS and the SuiteSparse library.
The simplest and recommended way is using the CMake package configuration
from the build directory (`AmiciConfig.cmake` for use with 
[`find_package`](https://cmake.org/cmake/help/latest/command/find_package.html)
which tells CMake about all AMICI dependencies.

During model import a [CMake](https://cmake.org/) file (`CMakeLists.txt`)
will be generated automatically. The CMake file shows the above-mentioned
library dependencies. `CMakeLists.txt`, together with `main.cpp`, provides a
scaffold for a standalone simulation program.
The required numerical libraries are shipped with AMICI.
To compile them, run `${AMICI_ROOT_DIR}/scripts/buildAll.sh` once. HDF5
libraries and header files need to be installed separately (see AMICI
installation instructions). 

More information on how to run the compiled sample program is provided in
`main.cpp`.

The file `wrapfunctions.h` exports some functions with generic names e.g. to
create a model instance. These functions are meant to be used by applications
which may be linked to single, potentially changing model to avoid hard-coding
the model name. Including multiple `wrapfunctions.h` files from different
models at once is not possible.
 

## Running a simulation

The complete AMICI API is available through `amici/amici.h`; this is the only
header file that needs to be included for basic usage.
Additionally, `amici/hdf5.h` may be helpful. This file provides some
functions for reading and writing [HDF5](https://support.hdfgroup.org/) files). 
All functions are declared within the `amici` namespace.

The entry function for running an AMICI simulation is
`amici::runAmiciSimulation(...)`, declared in `amici/amici.h`. 

This function requires 
    
* a `amici::Model` instance. For the example `model_steadystate` the respective
  class is provided as `Model_model_steadystate` in `model_steadystate.h`.
  For convenience, the header `wrapfunctions.h` defines a function
  `getModel()`, that returns an instance of that class.

* a `amici::Solver` instance. This solver instance needs to match the
  requirements of the model and can be generated using `model->getSolver()`.

* optionally an `amici::ExpData` instance, which contains any experimental data
  (e.g. measurements, noise model parameters or model inputs) to
  compute residuals or an objective function.

A scaffold for a standalone simulation program is generated in `main.cpp` in
the model source directory. This program shows how to use the
above-mentioned classes and how to obtain the simulation results.

For running simulations for multiple experimental conditions
(multiple `amici::ExpData` instances), `amici::runAmiciSimulations(...)`
provides an alternative entry point. If AMICI (and your application)
have been compiled with OpenMP support (see installation guide), this allows
for running those simulations in parallel.


## Parameter estimation for AMICI models in high-performance computing environments

To perform parameter estimation for large or otherwise computationally
demanding AMICI models from C++ in a high-performance computing environment,
you may find the [parPE](https://github.com/ICB-DCM/parPE/) library helpful.
parPE allows for the private or shared memory parallel evaluation of a cost
function requiring multiple simulations of the same model with different
inputs. It provides interfaces to different optimizers, such as Ipopt.
