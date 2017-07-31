#
# Build AMICI library
#

project(amici)

cmake_minimum_required(VERSION 2.8)

set(AMICI_DIR ${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")

add_definitions(-DAMICI_WITHOUT_MATLAB)

set(SUITESPARSE_DIR "${AMICI_DIR}/SuiteSparse/")
set(SUITESPARSE_LIB_DIR "${AMICI_DIR}/SuiteSparse/lib")
set(SUNDIALS_LIB_DIR "${AMICI_DIR}/sundials/build/lib")

find_package(HDF5 COMPONENTS C HL REQUIRED)

# AMICI requires BLAS, currently Intel MKL, CBLAS or MATLAB BLAS can be used.
# The latter is not supported via CMake yet.
set(BLAS "CBLAS" CACHE STRING "BLAS library to use")
set_property(CACHE BLAS PROPERTY STRINGS "CBLAS" "MKL")

if(${BLAS} STREQUAL "MKL")
    add_definitions(-DAMICI_BLAS_MKL)
    set(BLAS_LIBRARIES -lmkl CACHE STRING "")
    set(BLAS_INCLUDE_DIRS "" CACHE STRING "")
else()
    add_definitions(-DAMICI_BLAS_CBLAS)
    set(BLAS_INCLUDE_DIRS "" CACHE STRING "")
    set(BLAS_LIBRARIES -lcblas CACHE STRING "")
endif()

include_directories("${AMICI_DIR}")
include_directories("${HDF5_INCLUDE_DIRS}")
include_directories("${AMICI_DIR}/sundials/include")
include_directories("${SUITESPARSE_DIR}/include")

if(NOT "${BLAS_INCLUDE_DIRS}" STREQUAL "")
    include_directories(${BLAS_INCLUDE_DIRS})
endif()

# Library source files
set(AMICI_SRC_LIST
    ${AMICI_DIR}/src/symbolic_functions.cpp
    ${AMICI_DIR}/src/amici_interface_cpp.cpp
    ${AMICI_DIR}/src/amici_model_functions.cpp
    ${AMICI_DIR}/src/amici.cpp
    ${AMICI_DIR}/src/amici_misc.cpp
    ${AMICI_DIR}/src/udata.cpp
    ${AMICI_DIR}/src/rdata.cpp
    ${AMICI_DIR}/src/tdata.cpp
    ${AMICI_DIR}/src/edata.cpp
    ${AMICI_DIR}/src/amici_hdf5.cpp
    ${AMICI_DIR}/src/spline.cpp
)

add_library(${PROJECT_NAME} ${AMICI_SRC_LIST})

target_link_libraries(${PROJECT_NAME}
    ${SUNDIALS_LIB_DIR}/libsundials_nvecserial${SHARED_OBJECT_EXTENSION}
    ${SUNDIALS_LIB_DIR}/libsundials_cvodes${SHARED_OBJECT_EXTENSION}
    ${SUITESPARSE_LIB_DIR}/libcolamd${SHARED_OBJECT_EXTENSION}
    ${SUITESPARSE_LIB_DIR}/libklu${SHARED_OBJECT_EXTENSION}
    ${SUITESPARSE_LIB_DIR}/libbtf${SHARED_OBJECT_EXTENSION}
    ${SUITESPARSE_LIB_DIR}/libamd${SHARED_OBJECT_EXTENSION}
    ${SUITESPARSE_LIB_DIR}/libsuitesparseconfig${SHARED_OBJECT_EXTENSION}
    ${HDF5_HL_LIBRARIES}
    ${HDF5_C_LIBRARIES}
    ${BLAS_LIBRARIES}
    -ldl -lz -lm
)

