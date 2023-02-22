#
# Build AMICI library
#
cmake_minimum_required(VERSION 3.15)

project(amici)

message("CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")
message("CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
message("CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 14)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # require at least gcc 4.9, otherwise regex wont work properly
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9)
    message(FATAL_ERROR "GCC version must be at least 4.9!")
  endif()
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compiler flags
include(CheckCXXCompilerFlag)
set(MY_CXX_FLAGS -Wall)
foreach(FLAG ${MY_CXX_FLAGS})
  unset(CUR_FLAG_SUPPORTED CACHE)
  check_cxx_compiler_flag(${FLAG} CUR_FLAG_SUPPORTED)
  if(${CUR_FLAG_SUPPORTED})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAG}")
  endif()
endforeach(FLAG)

# Debug build?
if($ENV{ENABLE_AMICI_DEBUGGING} OR $ENV{ENABLE_GCOV_COVERAGE})
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -UNDEBUG")
    set(CMAKE_BUILD_TYPE "Debug")

    set(MY_CXX_FLAGS -Werror -Wno-error=deprecated-declarations)
    foreach(FLAG ${MY_CXX_FLAGS})
      unset(CUR_FLAG_SUPPORTED CACHE)
      check_cxx_compiler_flag(${FLAG} CUR_FLAG_SUPPORTED)
      if(${CUR_FLAG_SUPPORTED})
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAG}")
      endif()
    endforeach(FLAG)
endif()

# coverage options
if($ENV{ENABLE_GCOV_COVERAGE})
  string(APPEND CMAKE_CXX_FLAGS_DEBUG "--coverage")
  string(APPEND CMAKE_EXE_LINKER_FLAGS_DEBUG " --coverage")
endif()

if($ENV{AMICI_CXXFLAGS})
  string(APPEND CMAKE_CXX_FLAGS "$ENV{AMICI_CXXFLAGS}")
endif()
if($ENV{AMICI_LDFLAGS})
  string(APPEND CMAKE_EXE_LINKER_FLAGS "$ENV{AMICI_LDFLAGS}")
endif()

# misc options
option(AMICI_PYTHON_EXT_ONLY "Build only the Python extension?" OFF)

# find dependencies
include(GNUInstallDirs)

find_package(OpenMP)

option(AMICI_TRY_ENABLE_HDF5 "Build with HDF5 support if available?" ON)
option(ENABLE_HDF5 "Build with HDF5 support?" OFF)

if(ENABLE_HDF5)
  find_package(
    HDF5
    COMPONENTS C HL CXX
    REQUIRED)
elseif(AMICI_TRY_ENABLE_HDF5)
  find_package(HDF5 COMPONENTS C HL CXX)
endif()

set(VENDORED_SUNDIALS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/sundials)
set(VENDORED_SUNDIALS_BUILD_DIR ${VENDORED_SUNDIALS_DIR}/build)
set(VENDORED_SUNDIALS_INSTALL_DIR ${VENDORED_SUNDIALS_BUILD_DIR})
set(SUNDIALS_PRIVATE_INCLUDE_DIRS "${VENDORED_SUNDIALS_DIR}/src")
find_package(SUNDIALS REQUIRED PATHS
             "${VENDORED_SUNDIALS_INSTALL_DIR}/lib/cmake/sundials/")

set(GSL_LITE_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/gsl")

# AMICI requires BLAS, currently Intel MKL, CBLAS or MATLAB BLAS can be used.
# The latter is not supported via CMake yet.
set(BLAS
    "CBLAS"
    CACHE STRING "BLAS library to use")
set_property(CACHE BLAS PROPERTY STRINGS "CBLAS" "MKL" "ACCELERATE")
if(${BLAS} STREQUAL "MKL" OR DEFINED ENV{MKLROOT})
  if(DEFINED ENV{MKLROOT})
    # This is set by Environment Modules
    message(STATUS "Using MKL_INCDIR and MKL_LIB from environment module")
    set(BLAS
        "MKL"
        CACHE STRING "BLAS library to use" FORCE)
    set(BLAS_INCLUDE_DIRS
        "$ENV{MKL_INCDIR}"
        CACHE STRING "" FORCE)
    set(BLAS_LIBRARIES
        "$ENV{MKL_LIB}"
        CACHE STRING "" FORCE)
  else()
    set(BLAS_INCLUDE_DIRS
        ""
        CACHE STRING "")
    set(BLAS_LIBRARIES
        -lmkl
        CACHE STRING "")
  endif()
else()
  set(BLAS_INCLUDE_DIRS
      ""
      CACHE STRING "")
  set(BLAS_LIBRARIES
      -lcblas
      CACHE STRING "")
endif()
add_definitions(-DAMICI_BLAS_${BLAS})

# Add target to create version file
add_custom_target(
  version
  ${CMAKE_COMMAND}
  -D
  SRC=${CMAKE_CURRENT_SOURCE_DIR}/include/amici/version.in.h
  -D
  DST=${CMAKE_CURRENT_BINARY_DIR}/include/amici/version.h
  -P
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/configureVersion.cmake)

# Library source files
set(AMICI_SRC_LIST
    src/symbolic_functions.cpp
    src/cblas.cpp
    src/amici.cpp
    src/misc.cpp
    src/rdata.cpp
    src/edata.cpp
    src/exception.cpp
    src/simulation_parameters.cpp
    src/spline.cpp
    src/solver.cpp
    src/solver_cvodes.cpp
    src/solver_idas.cpp
    src/logging.cpp
    src/model.cpp
    src/model_ode.cpp
    src/model_dae.cpp
    src/model_state.cpp
    src/newton_solver.cpp
    src/forwardproblem.cpp
    src/steadystateproblem.cpp
    src/backwardproblem.cpp
    src/sundials_matrix_wrapper.cpp
    src/sundials_linsol_wrapper.cpp
    src/abstract_model.cpp
    src/vector.cpp
    include/amici/abstract_model.h
    include/amici/amici.h
    include/amici/backwardproblem.h
    include/amici/cblas.h
    include/amici/defines.h
    include/amici/edata.h
    include/amici/exception.h
    include/amici/forwardproblem.h
    include/amici/hdf5.h
    include/amici/logging.h
    include/amici/misc.h
    include/amici/model_dae.h
    include/amici/model_dimensions.h
    include/amici/model.h
    include/amici/model_ode.h
    include/amici/model_state.h
    include/amici/newton_solver.h
    include/amici/rdata.h
    include/amici/serialization.h
    include/amici/simulation_parameters.h
    include/amici/solver_cvodes.h
    include/amici/solver.h
    include/amici/solver_idas.h
    include/amici/spline.h
    include/amici/steadystateproblem.h
    include/amici/sundials_linsol_wrapper.h
    include/amici/sundials_matrix_wrapper.h
    include/amici/symbolic_functions.h
    include/amici/vector.h
    $<$<BOOL:${HDF5_FOUND}>:src/hdf5.cpp>)

add_library(${PROJECT_NAME} ${AMICI_SRC_LIST})

set(AMICI_CXX_OPTIONS
    ""
    CACHE STRING "C++ options for libamici (semicolon-separated)")
target_compile_options(${PROJECT_NAME} PRIVATE "${AMICI_CXX_OPTIONS}")

add_dependencies(${PROJECT_NAME} version)

file(
  GLOB PUBLIC_HEADERS
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  include/amici/*.h)
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER
                                                 "${PUBLIC_HEADERS}")
target_include_directories(
  ${PROJECT_NAME}
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
  PUBLIC $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PUBLIC $<BUILD_INTERFACE:${GSL_LITE_INCLUDE_DIR}>
  PUBLIC $<INSTALL_INTERFACE:${CMAKE_INSTALL_DATADIR}/amici/swig>
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/swig>
  PRIVATE ${SUNDIALS_PRIVATE_INCLUDE_DIRS})

if(NOT "${BLAS_INCLUDE_DIRS}" STREQUAL "")
  target_include_directories(${PROJECT_NAME} PUBLIC ${BLAS_INCLUDE_DIRS})
endif()

target_link_libraries(
  ${PROJECT_NAME}
  PUBLIC SUNDIALS::generic_static
  PUBLIC SUNDIALS::nvecserial_static
  PUBLIC SUNDIALS::sunmatrixband_static
  PUBLIC SUNDIALS::sunmatrixdense_static
  PUBLIC SUNDIALS::sunmatrixsparse_static
  PUBLIC SUNDIALS::sunlinsolband_static
  PUBLIC SUNDIALS::sunlinsoldense_static
  PUBLIC SUNDIALS::sunlinsolpcg_static
  PUBLIC SUNDIALS::sunlinsolspbcgs_static
  PUBLIC SUNDIALS::sunlinsolspfgmr_static
  PUBLIC SUNDIALS::sunlinsolspgmr_static
  PUBLIC SUNDIALS::sunlinsolsptfqmr_static
  PUBLIC SUNDIALS::sunlinsolklu_static
  PUBLIC SUNDIALS::sunnonlinsolnewton_static
  PUBLIC SUNDIALS::sunnonlinsolfixedpoint_static
  PUBLIC SUNDIALS::cvodes_static
  PUBLIC SUNDIALS::idas_static
  PUBLIC ${BLAS_LIBRARIES}
  PUBLIC ${CMAKE_DL_LIBS})

if(HDF5_FOUND)
  message(STATUS "HDF5 library found. Building AMICI with HDF5 support.")
  target_link_libraries(
    ${PROJECT_NAME}
    PUBLIC hdf5::hdf5_hl_cpp
    PUBLIC hdf5::hdf5_hl
    PUBLIC hdf5::hdf5_cpp
    PUBLIC hdf5::hdf5)
else()
  message(STATUS "HDF5 library NOT found. Building AMICI WITHOUT HDF5 support.")

endif()

if(OpenMP_FOUND)
  target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

option(SUNDIALS_SUPERLUMT_ENABLE "Enable sundials SuperLUMT?" OFF)
if(SUNDIALS_SUPERLUMT_ENABLE)
  target_link_libraries(${PROJECT_NAME}
                        PRIVATE SUNDIALS::sundials_sunlinsolsuperlumt)
endif()

# Create targets to make the sources show up in IDEs for convenience

# For matlab interface
if(NOT AMICI_PYTHON_EXT_ONLY)

  set(matlab_sources
      src/interface_matlab.cpp src/returndata_matlab.cpp
      include/amici/interface_matlab.h include/amici/returndata_matlab.h)
  find_package(Matlab)
  # In case we can find Matlab, we create a respective library to compile the
  # extension from cmake. Otherwise we just create a dummy target for the files
  # to show up inside IDEs. (Set the Matlab_ROOT_DIR cmake variable if CMake
  # cannot find your Matlab installation)
  if(${Matlab_FOUND})
    add_library(matlabInterface ${matlab_sources})
    set_target_properties(matlabInterface PROPERTIES INCLUDE_DIRECTORIES
                                                     "${Matlab_INCLUDE_DIRS}")

    target_link_libraries(matlabInterface PUBLIC amici)

  else()
    add_custom_target(matlabInterface SOURCES ${matlab_sources})
  endif()
  set_property(
    TARGET matlabInterface
    APPEND
    PROPERTY INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/include/")
endif()

# For template files
add_custom_target(
  fileTemplates
  SOURCES src/CMakeLists.template.cmake
          src/main.template.cpp
          src/model_header.ODE_template.h
          src/model.ODE_template.cpp
          src/wrapfunctions.ODE_template.h
          src/wrapfunctions.template.cpp
          swig/CMakeLists_model.cmake
          swig/modelname.template.i)
set_target_properties(
  fileTemplates PROPERTIES INCLUDE_DIRECTORIES
                           "${CMAKE_CURRENT_SOURCE_DIR}/include/")

if(NOT AMICI_PYTHON_EXT_ONLY)
  include(clang-tools)
  include(cmakelang-tools)
endif()

set(AUTHORS "Fabian Froehlich, Jan Hasenauer, Daniel Weindl and Paul Stapor")
set(AUTHOR_EMAIL "Fabian_Froehlich@hms.harvard.edu")

# <Export cmake configuration>
install(
  TARGETS ${PROJECT_NAME}
  EXPORT AmiciTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/amici)
export(
  EXPORT AmiciTargets
  FILE AmiciTargets.cmake
  NAMESPACE Upstream::)
include(CMakePackageConfigHelpers)
include(version)
configure_package_config_file(
  cmake/AmiciConfig.cmake "${CMAKE_CURRENT_BINARY_DIR}/AmiciConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Amici")
write_basic_package_version_file(AmiciConfigVersion.cmake
                                 COMPATIBILITY ExactVersion)
install(
  EXPORT AmiciTargets
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Amici"
  NAMESPACE Upstream::)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/AmiciConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/AmiciConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Amici)

install(DIRECTORY ThirdParty/gsl/gsl TYPE INCLUDE)
# When running from setup.py, this is a symlink we need to dereference
get_filename_component(_swig_realpath "swig" REALPATH)
install(DIRECTORY "${_swig_realpath}"
        DESTINATION ${CMAKE_INSTALL_DATADIR}/amici)

# Register package?
option(EXPORT_PACKAGE "Export AMICI library to CMake package registry?" ON)
if(EXPORT_PACKAGE)
  export(PACKAGE Amici)
endif()
# </Export cmake configuration>

# build interfaces for other languages
option(ENABLE_SWIG "Build AMICI swig library?" ON)
if(ENABLE_SWIG)
  add_subdirectory(swig)
endif()

option(ENABLE_PYTHON "Create Python module?" ON)
if(ENABLE_PYTHON AND NOT AMICI_PYTHON_EXT_ONLY)
  add_subdirectory(python)
endif()

option(BUILD_TESTS "Build integration tests?" ON)
if(BUILD_TESTS AND NOT AMICI_PYTHON_EXT_ONLY)
  if(HDF5_FOUND)
    enable_testing()

    add_subdirectory(tests/cpp)
  else()
    message(WARNING "Cannot build tests without HDF5 support.")
  endif()

endif()
