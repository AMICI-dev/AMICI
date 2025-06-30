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

if(DEFINED ENV{AMICI_BLAS_USE_SCIPY_OPENBLAS})
  message(
    STATUS
      "Using AMICI_BLAS_USE_SCIPY_OPENBLAS=${AMICI_BLAS_USE_SCIPY_OPENBLAS} from environment variable."
  )
  set(AMICI_BLAS_USE_SCIPY_OPENBLAS $ENV{AMICI_BLAS_USE_SCIPY_OPENBLAS})
endif()

if((${BLAS} STREQUAL "MKL" OR DEFINED ENV{MKLROOT})
   AND NOT AMICI_BLAS_USE_SCIPY_OPENBLAS)
  if(DEFINED ENV{MKLROOT})
    set(BLAS
        "MKL"
        CACHE STRING "BLAS library to use" FORCE)

    # Was MKLROOT set by /opt/intel/oneapi/setvars.sh? then cmake will find it
    message(STATUS "Trying to find BLAS based on MKLROOT=$ENV{MKLROOT}")
    # give the user the option to override the BLA_VENDOR
    if(NOT DEFINED BLA_VENDOR)
      set(BLA_VENDOR Intel10_64lp)
    endif()
    message(STATUS "Trying FindBLAS with BLA_VENDOR=${BLA_VENDOR}")
    find_package(BLAS)
    if(BLAS_FOUND)
      message(STATUS "Found BLAS via FindBLAS and MKLROOT")
    else()
      # This is set by Environment Modules and might not be compatible with
      # FindBLAS
      message(STATUS "Using MKL_INCDIR and MKL_LIB from environment module")
      set(BLAS_INCLUDE_DIRS
          "$ENV{MKL_INCDIR}"
          CACHE STRING "" FORCE)
      set(BLAS_LIBRARIES
          "$ENV{MKL_LIB}"
          CACHE STRING "" FORCE)
    endif()
  else()
    message(STATUS "BLAS is set to MKL, but MKLROOT is not set.")
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

  if(APPLE AND ((NOT DEFINED BLA_VENDOR OR BLA_VENDOR STREQUAL "All")
                AND NOT AMICI_BLAS_USE_SCIPY_OPENBLAS))
    set(BLA_VENDOR Apple)
    message(STATUS "Trying FindBLAS with BLA_VENDOR=${BLA_VENDOR}")
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
      set(AMICI_BLAS_USE_SCIPY_OPENBLAS FALSE)
    else()
      set(BLA_VENDOR "All")
    endif()
  endif()

  # Try the scipy-openblas64 package, assuming CMAKE_PREFIX_PATH is set to the
  # directory containing the package configuration.
  #
  # Set AMICI_BLAS_USE_SCIPY_OPENBLAS to TRUE if the package is found and used.
  # It must only be used in AmiciConfig.cmake if it was used for building AMICI.
  # CMAKE_PREFIX_PATH may contain the scipy-openblas64 directory, when building
  # model extensions, even if AMICI was not built with it originally.
  if(AMICI_BLAS_USE_SCIPY_OPENBLAS
     OR (NOT DEFINED AMICI_BLAS_USE_SCIPY_OPENBLAS
         AND (NOT DEFINED BLA_VENDOR OR BLA_VENDOR STREQUAL "All")))
    message(STATUS "Trying to find OpenBLAS in CONFIG mode (scipy-openblas64)")
    find_package(OpenBLAS CONFIG)
    if(OpenBLAS_FOUND)
      message(
        STATUS "Found OpenBLAS in CONFIG mode (OpenBLAS_DIR=${OpenBLAS_DIR})")
      set(BLAS_INCLUDE_DIRS ${OpenBLAS_INCLUDE_DIRS})
      set(BLAS_LIBRARIES ${OpenBLAS_LIBRARIES})
      # fix incorrect path, replace /bin/ by /lib/
      # https://github.com/MacPython/openblas-libs/issues/202
      string(REPLACE "/bin/" "/lib/" BLAS_LIBRARIES "${BLAS_LIBRARIES}")
      string(REPLACE "/libscipy_openblas64_.dll" "/libscipy_openblas64_.lib"
                     BLAS_LIBRARIES "${BLAS_LIBRARIES}")

      list(APPEND BLAS_DEFINES "BLAS_PREFIX=scipy_cblas_" "BLAS_SUFFIX=64_")
      set(BLAS_FOUND TRUE)
      set(AMICI_BLAS_USE_SCIPY_OPENBLAS TRUE)
    else()
      set(AMICI_BLAS_USE_SCIPY_OPENBLAS FALSE)
      message(STATUS "Could not find scipy-openblas64.")
    endif()
  else()
    set(AMICI_BLAS_USE_SCIPY_OPENBLAS FALSE)
    message(
      STATUS "Not looking for scipy-openblas64 because BLA_VENDOR=${BLA_VENDOR}"
    )
  endif()

  if(NOT BLAS_FOUND)
    message(STATUS "Trying FindBLAS with BLA_VENDOR=${BLA_VENDOR}")
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
  set_target_properties(
    BLAS
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${BLAS_INCLUDE_DIRS}"
               INTERFACE_LINK_LIBRARIES "${BLAS_LIBRARIES}"
               INTERFACE_COMPILE_DEFINITIONS "${BLAS_DEFINES}")
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
