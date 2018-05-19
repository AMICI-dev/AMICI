#
# Build AMICI library
#
cmake_minimum_required(VERSION 3.3)
project(amici)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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
find_package(HDF5 COMPONENTS C HL CXX REQUIRED)
set(HDF5_LIBRARIES ${HDF5_HL_LIBRARIES} ${HDF5_C_LIBRARIES} ${HDF5_CXX_LIBRARIES})
set(SUITESPARSE_DIR "${CMAKE_SOURCE_DIR}/ThirdParty/SuiteSparse/")
set(SUITESPARSE_INCLUDE_DIRS "${SUITESPARSE_DIR}/include" "${CMAKE_SOURCE_DIR}/ThirdParty/sundials/src")
set(SUITESPARSE_LIBRARIES
    ${SUITESPARSE_DIR}/KLU/Lib/libklu${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUITESPARSE_DIR}/COLAMD/Lib/libcolamd${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUITESPARSE_DIR}/BTF/Lib/libbtf${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUITESPARSE_DIR}/AMD/Lib/libamd${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUITESPARSE_DIR}/SuiteSparse_config/libsuitesparseconfig${CMAKE_STATIC_LIBRARY_SUFFIX}
    )

set(SUNDIALS_LIB_DIR "${CMAKE_SOURCE_DIR}/ThirdParty/sundials/build/lib")
set(SUNDIALS_LIBRARIES ${SUNDIALS_LIB_DIR}/libsundials_nvecserial${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUNDIALS_LIB_DIR}/libsundials_cvodes${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUNDIALS_LIB_DIR}/libsundials_idas${CMAKE_STATIC_LIBRARY_SUFFIX}
    )
set(SUNDIALS_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/ThirdParty/sundials/include")

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
    ${CMAKE_SOURCE_DIR}/src/hdf5.cpp
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
    )

add_library(${PROJECT_NAME} ${AMICI_SRC_LIST})
add_dependencies(${PROJECT_NAME} version)
file(GLOB PUBLIC_HEADERS include/amici/*.h)
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${PUBLIC_HEADERS}")
target_include_directories(${PROJECT_NAME}
    PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
    PUBLIC swig
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
add_custom_target(matlabInterface SOURCES src/interface_matlab.cpp src/returndata_matlab.cpp)
add_custom_target(fileTemplates SOURCES src/CMakeLists.template.txt src/main.template.cpp src/wrapfunctions.ODE_template.h src/wrapfunctions.template.cpp swig/CMakeLists_model.txt)

if($ENV{ENABLE_GCOV_COVERAGE})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 --coverage")
endif()

include(clang-tools)

set(AUTHORS "Fabian Froehlich, Jan Hasenauer, Daniel Weindl and Paul Stapor")
set(AUTHOR_EMAIL "fabian.froehlich@helmholtz-muenchen.de")

# <Export cmake configuration>
include(GNUInstallDirs)
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
get_property(AmiciConfigIncludes DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
list(APPEND AmiciConfigIncludes ${CMAKE_CURRENT_SOURCE_DIR}/swig/)
configure_file(cmake/AmiciConfig.cmake
    "${CMAKE_CURRENT_BINARY_DIR}/AmiciConfig.cmake"
    @ONLY
    )
include(CMakePackageConfigHelpers)
include(version)
write_basic_package_version_file(AmiciConfigVersion.cmake COMPATIBILITY ExactVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/AmiciConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/AmiciTargets.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/AmiciConfigVersion.cmake
    DESTINATION share/Amici/cmake )
# Register package
export(PACKAGE Amici)
# </Export cmake configuration>


# build interfaces for other languages
add_subdirectory(swig)
add_subdirectory(python)
