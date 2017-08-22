project(amiciIntegrationTests)

cmake_minimum_required(VERSION 2.8)

set(CMAKE_BUILD_TYPE Debug)

set(AMICI_DIR "${CMAKE_SOURCE_DIR}/../../")
set(AMICI_LIBRARIES "${AMICI_DIR}/build/libamici.a")
set(BLAS_LIBRARIES "-lcblas")
set(SUITESPARSE_DIR "${AMICI_DIR}/SuiteSparse/")
set(SUITESPARSE_LIB_DIR "${AMICI_DIR}/SuiteSparse/lib")
set(SUNDIALS_LIB_DIR "${AMICI_DIR}/sundials/build/lib")
set(CPPUTEST_DIR "${AMICI_DIR}/ThirdParty/cpputest-3.8/")
set(CPPUTEST_INCLUDE_DIR "${CPPUTEST_DIR}/include")
set(CPPUTEST_LIBRARY "${CPPUTEST_DIR}/cpputest_build/lib/libCppUTest.a" "${CPPUTEST_DIR}/cpputest_build/lib/libCppUTestExt.a")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -include ${CPPUTEST_INCLUDE_DIR}/CppUTest/MemoryLeakDetectorNewMacros.h")

find_package(HDF5 COMPONENTS C HL REQUIRED)

include_directories("${AMICI_DIR}")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories("${HDF5_INCLUDE_DIRS}")
include_directories("${AMICI_DIR}/sundials/include")
include_directories("${SUITESPARSE_DIR}/include")
include_directories(${CPPUTEST_INCLUDE_DIR})

if(APPLE)
    set(SHARED_OBJECT_EXTENSION .dylib)
else()
    set(SHARED_OBJECT_EXTENSION .so)
endif()

add_library(amici STATIC IMPORTED)
set_property(TARGET amici PROPERTY IMPORTED_LOCATION "${AMICI_DIR}/build/libamici.a")
set_property(TARGET amici PROPERTY INTERFACE_LINK_LIBRARIES
    "${SUNDIALS_LIB_DIR}/libsundials_nvecserial${SHARED_OBJECT_EXTENSION}"
    "${SUNDIALS_LIB_DIR}/libsundials_cvodes${SHARED_OBJECT_EXTENSION}"
    "${SUITESPARSE_LIB_DIR}/libcolamd${SHARED_OBJECT_EXTENSION}"
    "${SUITESPARSE_LIB_DIR}/libklu${SHARED_OBJECT_EXTENSION}"
    "${SUITESPARSE_LIB_DIR}/libbtf${SHARED_OBJECT_EXTENSION}"
    "${SUITESPARSE_LIB_DIR}/libamd${SHARED_OBJECT_EXTENSION}"
    "${SUITESPARSE_LIB_DIR}/libsuitesparseconfig${SHARED_OBJECT_EXTENSION}"
    "${HDF5_HL_LIBRARIES}"
    "${HDF5_C_LIBRARIES}"
    "${BLAS_LIBRARIES}"
    "-lpthread -ldl -lz"
    "-lm"
    )


enable_testing()

add_subdirectory(dirac)
add_subdirectory(steadystate)
add_subdirectory(jakstat_adjoint)
add_subdirectory(neuron)
add_subdirectory(jakstat_adjoint_o2)
add_subdirectory(neuron_o2)
