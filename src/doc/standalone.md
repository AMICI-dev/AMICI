# Using AMICI-generated code outside Matlab

AMICI (amiwrap.m)  translates the model definition into C++ code which is then compiled into a mex file for MATLAB. Advanced users can use this code within stand-online C/C++ applications for use in non-MATLAB environments (e.g. on high performance computing systems). This section will give a short overview over the generated files and provide a brief introduction of how this code can be included in other applications.


## Generated model files
amiwrap.m usually write the model source files to ${AMICI_ROOT_DIR}/models/${MODEL_NAME} by default. 
The content of a model source directory might look something like this (given `MODEL_NAME=model_steadystate`): 

```
CMakeLists.txt
hashes.mat 
main.cpp 
model_steadystate_deltaqB.cpp 
model_steadystate_deltaqB.h 
model_steadystate_deltaqB_mexa64.md5 
model_steadystate_deltaqB.o 
[... many more files model_steadystate_*.(cpp|h|md5|o) ]
wrapfunctions.cpp 
wrapfunctions.h 
wrapfunctions.o
```

Only `*.cpp` and `*.h` files will be needed for the model; `*.o` and `*.md5` are not required. 

## Running a simulation

The entry function for running an AMICI simulation is getSimulationResults(), declared in amici.h. This function requires all AMICI options and any experimental data. All options that would normally be passed to `simulate_${MODEL_NAME}()` in MATLAB are passed in a UserData struct (see `udata.h` for info). Any experimental data will be passed as ExpData struct (`edata.h`). The simulation results will be returned in a ReturnData struct (see `rdata.h`).

A scaffold for a standalone simulation program is generated in `main.cpp` in the model source directory. This programm shows how to initialize the above-mentioned structs and how to obtain the simulation results.

## Compiling and linking

The complete AMICI API is available through amici.h; this is the only header file that needs to be included. (There are some accessor macro definitions available in udata_accessors.h, rdata_accessors.h and edata_accessors.h which provide shortcuts for accessing struct members of UserData, ReturnData, ExpData, respectively. `ami_hdf5.h` provides some functions for reading and writing [HDF5](https://support.hdfgroup.org/) files). 

You need to compile and link `${AMICI_ROOT_DIR}/models/${MODEL_NAME}/*.cpp`, `${AMICI_ROOT_DIR}/src/*.cpp`, the SUNDIALS and the SUITESPARSE library.

Along with `main.cpp`, a [CMake](https://cmake.org/) file (`CMakeLists.txt`) will be generated automatically. The CMake file shows the abovementioned library dependencies. These files provide a scaffold for a standalone simulation program. (NOTE: This program should compile and link, but will crash most certainly without further problem-specific adaptations.)
