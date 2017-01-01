# C++ Interface

The @ref python and matlab @ref matlab interfaces can translate the model definition into C++ code, which is then compiled into a .mex file or a python module. Advanced users can also use this code within stand-alone C/C++ application for use in other environments (e.g. on high performance computing systems). This section will give a short overview over the generated files and provide a brief introduction of how this code can be included in other applications.

## Generated model files
amiwrap.m usually write the model source files to ${AMICI_ROOT_DIR}/models/${MODEL_NAME} by default. 
The content of a model source directory might look something like this (given `MODEL_NAME=model_steadystate`): 

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

Only `*.cpp` and `*.h` files will be needed for the model; `*.o` and `*.md5` are not required. 

## Running a simulation

The entry function for running an AMICI simulation is `runAmiciSimulation(...)`, declared in amici.h. This function requires 
    (i) a `Model` instance. For the example `model_steadystate` the respective class is provided as `Model_model_steadystate`  in `model_steadystate.h`. For convenience, the header `wrapfunctions.h` defines a function `getModel()`, that returns an instance of that class.
    (ii) a `Solver` instance. This solver instance needs to match the requirements of the model and can be generated using `model->getSolver()`.
    (iii) optionally an `ExpData` instance, which contains any experimental data.

A scaffold for a standalone simulation program is generated in `main.cpp` in the model source directory. This programm shows how to initialize the above-mentioned structs and how to obtain the simulation results.

## Compiling and linking

The complete AMICI API is available through `amici.h`; this is the only header file that needs to be included.  `hdf5.h` provides some functions for reading and writing [HDF5](https://support.hdfgroup.org/) files). 

You need to compile and link `${AMICI_ROOT_DIR}/models/${MODEL_NAME}/*.cpp`,  `${AMICI_ROOT_DIR}/src/*.cpp`, the SUNDIALS and the SUITESPARSE library.

Along with `main.cpp`, a [CMake](https://cmake.org/) file (`CMakeLists.txt`) will be generated automatically. The CMake file shows the abovementioned library dependencies. These files provide a scaffold for a standalone simulation program. The required numerical libraries are shipped with AMICI. To compile them, run `${AMICI_ROOT_DIR}/scripts/run-tests.sh` once. HDF5 libraries and header files need to be installed separately. 
More information on how to run the compiled program is provided in `main.cpp`.
