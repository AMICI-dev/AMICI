#
# Build AMICI library
#
cmake_minimum_required(VERSION 3.3)

if(POLICY CMP0065)
  cmake_policy(SET CMP0065 NEW)
endif(POLICY CMP0065)
if(POLICY CMP0077)
  cmake_policy(SET CMP0077 NEW)
endif(POLICY CMP0077)
if(POLICY CMP0074)
  # Use package_ROOT environment variables
  cmake_policy(SET CMP0074 NEW)
endif(POLICY CMP0074)

project(amici)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 14)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # require at least gcc 4.9, otherwise regex wont work properly
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9)
        message(FATAL_ERROR "GCC version must be at least 4.9!")
    endif()
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(ENABLE_SWIG ON)

# Compiler flags
include(CheckCXXCompilerFlag)
set(MY_CXX_FLAGS -Wall)
foreach(FLAG ${MY_CXX_FLAGS})
    unset(CUR_FLAG_SUPPORTED CACHE)
    CHECK_CXX_COMPILER_FLAG(${FLAG} CUR_FLAG_SUPPORTED)
    if(${CUR_FLAG_SUPPORTED})
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAG}")
    endif()
endforeach(FLAG)

# find dependencies
include(GNUInstallDirs)

option(ENABLE_HDF5 "Build with HDF5 support?" ON)
if(ENABLE_HDF5)
    find_package(HDF5 COMPONENTS C HL CXX REQUIRED)
    set(HDF5_LIBRARIES ${HDF5_HL_LIBRARIES} ${HDF5_C_LIBRARIES} ${HDF5_CXX_LIBRARIES})
endif()

set(SUITESPARSE_DIR "${CMAKE_SOURCE_DIR}/ThirdParty/SuiteSparse/")
set(SUITESPARSE_INCLUDE_DIRS "${SUITESPARSE_DIR}/include" "${CMAKE_SOURCE_DIR}/ThirdParty/sundials/src")
set(SUITESPARSE_LIBRARIES
    ${SUITESPARSE_DIR}/KLU/Lib/libklu${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUITESPARSE_DIR}/COLAMD/Lib/libcolamd${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUITESPARSE_DIR}/BTF/Lib/libbtf${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUITESPARSE_DIR}/AMD/Lib/libamd${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUITESPARSE_DIR}/SuiteSparse_config/libsuitesparseconfig${CMAKE_STATIC_LIBRARY_SUFFIX}
    )

set(SUNDIALS_LIB_DIR "${CMAKE_SOURCE_DIR}/ThirdParty/sundials/build/${CMAKE_INSTALL_LIBDIR}")
set(SUNDIALS_LIBRARIES
    ${SUNDIALS_LIB_DIR}/libsundials_nvecserial${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUNDIALS_LIB_DIR}/libsundials_sunlinsolband${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUNDIALS_LIB_DIR}/libsundials_sunlinsolklu${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUNDIALS_LIB_DIR}/libsundials_sunlinsolpcg${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUNDIALS_LIB_DIR}/libsundials_sunlinsolspbcgs${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUNDIALS_LIB_DIR}/libsundials_sunlinsolspfgmr${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUNDIALS_LIB_DIR}/libsundials_sunmatrixband${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUNDIALS_LIB_DIR}/libsundials_sunmatrixdense${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUNDIALS_LIB_DIR}/libsundials_sunmatrixsparse${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUNDIALS_LIB_DIR}/libsundials_sunnonlinsolfixedpoint${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUNDIALS_LIB_DIR}/libsundials_sunnonlinsolnewton${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUNDIALS_LIB_DIR}/libsundials_cvodes${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUNDIALS_LIB_DIR}/libsundials_idas${CMAKE_STATIC_LIBRARY_SUFFIX}
    )
set(SUNDIALS_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/ThirdParty/sundials/build/include")

option(SUNDIALS_SUPERLUMT_ENABLE "Enable sundials SuperLUMT?" OFF)
if(SUNDIALS_SUPERLUMT_ENABLE)
    set(SUNDIALS_LIBRARIES ${SUNDIALS_LIBRARIES}
        ${SUNDIALS_LIB_DIR}/libsundials_sunlinsolsuperlumt${CMAKE_STATIC_LIBRARY_SUFFIX}
        ${CMAKE_SOURCE_DIR}/ThirdParty/SuperLU_MT_3.1/lib/libsuperlu_mt_PTHREAD${CMAKE_STATIC_LIBRARY_SUFFIX}
        -lblas
        )
    set(SUNDIALS_INCLUDE_DIRS ${SUNDIALS_INCLUDE_DIRS}
        "${CMAKE_SOURCE_DIR}/ThirdParty/SuperLU_MT_3.1/SRC/")
endif()

set(GSL_LITE_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/ThirdParty/gsl")

# AMICI requires BLAS, currently Intel MKL, CBLAS or MATLAB BLAS can be used.
# The latter is not supported via CMake yet.
set(BLAS "CBLAS" CACHE STRING "BLAS library to use")
set_property(CACHE BLAS PROPERTY STRINGS "CBLAS" "MKL" "ACCELERATE")
if(${BLAS} STREQUAL "MKL" OR DEFINED ENV{MKLROOT})
    if(DEFINED ENV{MKLROOT})
        # This is set by Environment Modules
        message(STATUS "Using MKL_INCDIR and MKL_LIB from environment module")
        set(BLAS "MKL" CACHE STRING "BLAS library to use" FORCE)
        set(BLAS_INCLUDE_DIRS "$ENV{MKL_INCDIR}" CACHE STRING "" FORCE)
        set(BLAS_LIBRARIES "$ENV{MKL_LIB}" CACHE STRING "" FORCE)
    else()
        set(BLAS_INCLUDE_DIRS "" CACHE STRING "")
        set(BLAS_LIBRARIES -lmkl CACHE STRING "")
    endif()
else()
    set(BLAS_INCLUDE_DIRS "" CACHE STRING "")
    set(BLAS_LIBRARIES -lcblas CACHE STRING "")
endif()
add_definitions(-DAMICI_BLAS_${BLAS})

# Add target to create version file
add_custom_target(
    version
    ${CMAKE_COMMAND} -D SRC=${CMAKE_SOURCE_DIR}/include/amici/version.in.h
    -D DST=${CMAKE_BINARY_DIR}/include/amici/version.h
    -P ${CMAKE_SOURCE_DIR}/cmake/configureVersion.cmake
    )

# Library source files
set(AMICI_SRC_LIST
    ${CMAKE_SOURCE_DIR}/src/symbolic_functions.cpp
    ${CMAKE_SOURCE_DIR}/src/cblas.cpp
    ${CMAKE_SOURCE_DIR}/src/amici.cpp
    ${CMAKE_SOURCE_DIR}/src/misc.cpp
    ${CMAKE_SOURCE_DIR}/src/rdata.cpp
    ${CMAKE_SOURCE_DIR}/src/edata.cpp
    ${CMAKE_SOURCE_DIR}/src/exception.cpp
    ${CMAKE_SOURCE_DIR}/src/simulation_parameters.cpp
    ${CMAKE_SOURCE_DIR}/src/spline.cpp
    ${CMAKE_SOURCE_DIR}/src/solver.cpp
    ${CMAKE_SOURCE_DIR}/src/solver_cvodes.cpp
    ${CMAKE_SOURCE_DIR}/src/solver_idas.cpp
    ${CMAKE_SOURCE_DIR}/src/model.cpp
    ${CMAKE_SOURCE_DIR}/src/model_ode.cpp
    ${CMAKE_SOURCE_DIR}/src/model_dae.cpp
    ${CMAKE_SOURCE_DIR}/src/newton_solver.cpp
    ${CMAKE_SOURCE_DIR}/src/forwardproblem.cpp
    ${CMAKE_SOURCE_DIR}/src/steadystateproblem.cpp
    ${CMAKE_SOURCE_DIR}/src/backwardproblem.cpp
    ${CMAKE_SOURCE_DIR}/src/sundials_matrix_wrapper.cpp
    ${CMAKE_SOURCE_DIR}/src/sundials_linsol_wrapper.cpp
    ${CMAKE_SOURCE_DIR}/src/abstract_model.cpp
    ${CMAKE_SOURCE_DIR}/src/vector.cpp
    )
if(ENABLE_HDF5)
    list(APPEND AMICI_SRC_LIST ${CMAKE_SOURCE_DIR}/src/hdf5.cpp)
endif()

add_library(${PROJECT_NAME} ${AMICI_SRC_LIST})
add_dependencies(${PROJECT_NAME} version)
file(GLOB PUBLIC_HEADERS include/amici/*.h)
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${PUBLIC_HEADERS}")
target_include_directories(${PROJECT_NAME}
    PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
    PUBLIC swig
    PUBLIC ${GSL_LITE_INCLUDE_DIR}
    PUBLIC ${SUNDIALS_INCLUDE_DIRS}
    PUBLIC ${SUITESPARSE_INCLUDE_DIRS}
    PUBLIC ${HDF5_INCLUDE_DIRS}
    )

if(NOT "${BLAS_INCLUDE_DIRS}" STREQUAL "")
    target_include_directories(${PROJECT_NAME} PUBLIC ${BLAS_INCLUDE_DIRS})
endif()

target_link_libraries(${PROJECT_NAME}
    PUBLIC ${SUNDIALS_LIBRARIES}
    PUBLIC ${SUITESPARSE_LIBRARIES}
    PUBLIC ${HDF5_LIBRARIES}
    PUBLIC ${BLAS_LIBRARIES}
    )

# Create targets to make the sources show up in IDEs for convenience

# For matlab interface
add_custom_target(matlabInterface
    SOURCES
    src/interface_matlab.cpp
    src/returndata_matlab.cpp
    include/amici/interface_matlab.h)
set_target_properties(matlabInterface
    PROPERTIES INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/include/")
find_package(Matlab)
if (${Matlab_FOUND})
    # In case we can find Matlab, use the respective include directories
    # for better IDE integration (set Matlab_ROOT_DIR cmake variable
    # if CMake cannot find your Matlab installation)
    set_property(TARGET matlabInterface APPEND
        PROPERTY INCLUDE_DIRECTORIES "${Matlab_INCLUDE_DIRS}")
endif()

# For template files
add_custom_target(
    fileTemplates
    SOURCES
    src/CMakeLists.template.cmake
    src/main.template.cpp
    src/model_header.ODE_template.h
    src/model.ODE_template.cpp
    src/wrapfunctions.ODE_template.h
    src/wrapfunctions.template.cpp
    swig/CMakeLists_model.cmake
    swig/modelname.template.i
    )
set_target_properties(fileTemplates
    PROPERTIES INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/include/")


if($ENV{ENABLE_GCOV_COVERAGE})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 --coverage")
endif()

include(clang-tools)

set(AUTHORS "Fabian Froehlich, Jan Hasenauer, Daniel Weindl and Paul Stapor")
set(AUTHOR_EMAIL "Fabian_Froehlich@hms.harvard.edu")

# <Export cmake configuration>
install(TARGETS ${PROJECT_NAME} EXPORT AmiciTargets
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/amici
)
export(EXPORT AmiciTargets FILE AmiciTargets.cmake
    NAMESPACE Upstream::
    )
include(CMakePackageConfigHelpers)
include(version)
configure_package_config_file(
    cmake/AmiciConfig.cmake
    "${CMAKE_CURRENT_BINARY_DIR}/AmiciConfig.cmake"
    INSTALL_DESTINATION "${LIB_INSTALL_DIR}/cmake/"
    )
write_basic_package_version_file(AmiciConfigVersion.cmake COMPATIBILITY ExactVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/AmiciConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/AmiciTargets.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/AmiciConfigVersion.cmake
    DESTINATION share/Amici/cmake )
# Register package

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
if(ENABLE_PYTHON)
    add_subdirectory(python)
endif()

option(BUILD_TESTS "Build integration tests?" ON)
if(BUILD_TESTS)
    if(ENABLE_HDF5)
        enable_testing()

        add_subdirectory(tests/cpputest)
    else()
        message(WARNING "Cannot build tests with ENABLE_HDF5=OFF.")
    endif()

endif()
