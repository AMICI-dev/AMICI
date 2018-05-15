#
# Build AMICI library
#

project(amici)

cmake_minimum_required(VERSION 2.8)

set(AMICI_DIR ${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(CheckCXXCompilerFlag)
set(MY_CXX_FLAGS -Wall)
foreach(FLAG ${MY_CXX_FLAGS})
    unset(CUR_FLAG_SUPPORTED CACHE)
    CHECK_CXX_COMPILER_FLAG(${FLAG} CUR_FLAG_SUPPORTED)
    if(${CUR_FLAG_SUPPORTED})
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAG}")
    endif()
endforeach(FLAG)

add_definitions(-DAMICI_WITHOUT_MATLAB)

set(SUITESPARSE_DIR "${AMICI_DIR}/SuiteSparse/")
set(SUNDIALS_LIB_DIR "${AMICI_DIR}/sundials/build/lib")
set(CMAKE_MODULE_PATH ${AMICI_DIR}/CMakeModules)

find_package(HDF5 COMPONENTS C HL CXX REQUIRED)

# AMICI requires BLAS, currently Intel MKL, CBLAS or MATLAB BLAS can be used.
# The latter is not supported via CMake yet.
set(BLAS "CBLAS" CACHE STRING "BLAS library to use")
set_property(CACHE BLAS PROPERTY STRINGS "CBLAS" "MKL" "ACCELERATE")
add_definitions(-DAMICI_BLAS_${BLAS})
if(${BLAS} STREQUAL "MKL")
    set(BLAS_INCLUDE_DIRS "" CACHE STRING "")
    set(BLAS_LIBRARIES -lmkl CACHE STRING "")
else()
    set(BLAS_INCLUDE_DIRS "" CACHE STRING "")
    set(BLAS_LIBRARIES -lcblas CACHE STRING "")
endif()

include_directories("${AMICI_DIR}/include/")
include_directories("${HDF5_INCLUDE_DIRS}")
include_directories("${AMICI_DIR}/sundials/include")
include_directories("${AMICI_DIR}/sundials/src")
include_directories("${SUITESPARSE_DIR}/include")

if(NOT "${BLAS_INCLUDE_DIRS}" STREQUAL "")
    include_directories(${BLAS_INCLUDE_DIRS})
endif()

# Library source files
set(AMICI_SRC_LIST
    ${AMICI_DIR}/src/symbolic_functions.cpp
    ${AMICI_DIR}/src/cblas.cpp
    ${AMICI_DIR}/src/amici.cpp
    ${AMICI_DIR}/src/misc.cpp
    ${AMICI_DIR}/src/rdata.cpp
    ${AMICI_DIR}/src/edata.cpp
    ${AMICI_DIR}/src/hdf5.cpp
    ${AMICI_DIR}/src/spline.cpp
    ${AMICI_DIR}/src/solver.cpp
    ${AMICI_DIR}/src/solver_cvodes.cpp
    ${AMICI_DIR}/src/solver_idas.cpp
    ${AMICI_DIR}/src/model.cpp
    ${AMICI_DIR}/src/model_ode.cpp
    ${AMICI_DIR}/src/model_dae.cpp
    ${AMICI_DIR}/src/newton_solver.cpp
    ${AMICI_DIR}/src/forwardproblem.cpp
    ${AMICI_DIR}/src/steadystateproblem.cpp
    ${AMICI_DIR}/src/backwardproblem.cpp
)

add_library(${PROJECT_NAME} ${AMICI_SRC_LIST})

target_link_libraries(${PROJECT_NAME}
    ${SUNDIALS_LIB_DIR}/libsundials_nvecserial${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUNDIALS_LIB_DIR}/libsundials_cvodes${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUNDIALS_LIB_DIR}/libsundials_idas${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUITESPARSE_DIR}/KLU/Lib/libklu${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUITESPARSE_DIR}/COLAMD/Lib/libcolamd${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUITESPARSE_DIR}/BTF/Lib/libbtf${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUITESPARSE_DIR}/AMD/Lib/libamd${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUITESPARSE_DIR}/SuiteSparse_config/libsuitesparseconfig${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${HDF5_HL_LIBRARIES}
    ${HDF5_C_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
    ${BLAS_LIBRARIES}
    -ldl -lz -lm
)

# Create targets to make the sources show up in IDEs for convenience
add_custom_target(matlabInterface SOURCES src/interface_matlab.cpp src/returndata_matlab.cpp)
add_custom_target(fileTemplates SOURCES src/CMakeLists.template.txt src/main.template.cpp)

if($ENV{ENABLE_GCOV_COVERAGE})
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 --coverage")
endif()

include(CMakeModules/clang-tools.cmake)

# <Export cmake configuration>
include(GenerateExportHeader)
generate_export_header(amici)
install(TARGETS ${PROJECT_NAME} EXPORT AmiciTargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include
)
export(EXPORT AmiciTargets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/AmiciTargets.cmake"
  NAMESPACE Upstream::
)
get_property(AmiciConfigIncludes DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
configure_file(CMakeModules/AmiciConfig.cmake
  "${CMAKE_CURRENT_BINARY_DIR}/AmiciConfig.cmake"
  @ONLY
)
# Register package
export(PACKAGE Amici)
# </Export cmake configuration>
