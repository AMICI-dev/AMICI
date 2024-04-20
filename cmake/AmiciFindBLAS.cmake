# Find a BLAS
#
# AMICI requires BLAS, currently Intel MKL, CBLAS or MATLAB BLAS can be used.
# The latter is not supported via CMake yet.
#
# This module defines the BLAS::BLAS IMPORTED target.

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
elseif(
  NOT DEFINED ENV{BLAS_LIBS}
  AND NOT DEFINED ENV{BLAS_CFLAGS}
  AND NOT BLAS_FOUND)
  # if nothing is specified via environment variables, let's try FindBLAS
  if($(CMAKE_VERSION) VERSION_GREATER_EQUAL 3.22)
    set(BLA_SIZEOF_INTEGER 8)
  endif()

  if(APPLE AND (NOT DEFINED BLA_VENDOR OR BLA_VENDOR STREQUAL "All"))
    set(BLA_VENDOR Apple)
    find_package(BLAS)
    if(BLAS_FOUND)
      message(STATUS "Found Apple Accelerate BLAS")
      set_property(
        TARGET BLAS::BLAS
        APPEND
        PROPERTY INTERFACE_COMPILE_DEFINITIONS ACCELERATE_NEW_LAPACK)
      set_property(
        TARGET BLAS::BLAS
        APPEND
        PROPERTY INTERFACE_COMPILE_DEFINITIONS ACCELERATE_LAPACK_ILP64)
    else()
      set(BLA_VENDOR "All")
    endif()

  endif()
  if(NOT BLAS_FOUND)
    find_package(BLAS)
    if(BLAS_FOUND)
      message(STATUS "Found BLAS via FindBLAS")
    endif()
  endif()
  if(NOT BLAS_FOUND)
    # Nothing specified by the user and FindBLAS didn't find anything; let's try
    # if cblas is available on the system paths.
    set(BLAS_INCLUDE_DIRS
        ""
        CACHE STRING "")
    set(BLAS_LIBRARIES
        -lcblas
        CACHE STRING "")
  endif()
endif()

# Create an imported target
if(NOT TARGET BLAS::BLAS)
  add_library(BLAS INTERFACE)
  set_target_properties(BLAS PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${BLAS_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${BLAS_LIBRARIES}")
  add_library(BLAS::BLAS ALIAS BLAS)
  if("${PROJECT_NAME}" STREQUAL "amici")
      install(TARGETS BLAS EXPORT BLAS)
      export(EXPORT BLAS NAMESPACE BLAS::)
      install(
        EXPORT BLAS
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Amici"
        NAMESPACE BLAS::)
  endif()

  # legacy python package environment variables:
  if(DEFINED ENV{BLAS_CFLAGS})
    target_compile_options(BLAS INTERFACE "$ENV{BLAS_CFLAGS}")
  endif()
  if(DEFINED ENV{BLAS_LIBS})
    # Note that, on Windows, at least for ninja, this will only work with dashes
    # instead of slashes in any linker options
    target_link_libraries(BLAS INTERFACE "$ENV{BLAS_LIBS}")
  endif()

  if(NOT "${BLAS_INCLUDE_DIRS}" STREQUAL "")
    target_include_directories(BLAS INTERFACE ${BLAS_INCLUDE_DIRS})
  endif()
  target_compile_definitions(BLAS INTERFACE "AMICI_BLAS_${BLAS}")
else()
  target_compile_definitions(BLAS::BLAS INTERFACE "AMICI_BLAS_${BLAS}")
endif()
