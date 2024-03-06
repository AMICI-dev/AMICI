@PACKAGE_INIT@

include(CMakeFindDependencyMacro)

find_package(OpenMP)

# Current SUNDIALSConfig.cmake doesn't use KLU's FindKLU, but hardcodes paths,
# and is not relocatable. This does not work with Python package installation in
# tmpdirs.
list(APPEND CMAKE_MODULE_PATH
     @CMAKE_SOURCE_DIR@/ThirdParty/SuiteSparse/lib/cmake/SuiteSparse/)
if(NOT DEFINED SUITESPARSE_CONFIG_ROOT)
  set(SUITESPARSE_CONFIG_ROOT @CMAKE_SOURCE_DIR@/ThirdParty/SuiteSparse/)
endif()
if(NOT DEFINED AMD_ROOT)
  set(AMD_ROOT @CMAKE_SOURCE_DIR@/ThirdParty/SuiteSparse/)
endif()
if(NOT DEFINED BTF_ROOT)
  set(BTF_ROOT @CMAKE_SOURCE_DIR@/ThirdParty/SuiteSparse/)
endif()
if(NOT DEFINED COLAMD_ROOT)
  set(COLAMD_ROOT @CMAKE_SOURCE_DIR@/ThirdParty/SuiteSparse/)
endif()
if(NOT DEFINED KLU_ROOT)
  set(KLU_ROOT @CMAKE_SOURCE_DIR@/ThirdParty/SuiteSparse/)
endif()
find_dependency(SuiteSparse_config REQUIRED)
find_dependency(AMD REQUIRED)
find_dependency(BTF REQUIRED)
find_dependency(COLAMD REQUIRED)
find_dependency(KLU REQUIRED)

include("${CMAKE_CURRENT_LIST_DIR}/AmiciFindBLAS.cmake")

add_library(SUNDIALS::KLU INTERFACE IMPORTED)
target_link_libraries(
  SUNDIALS::KLU
  INTERFACE "${KLU_STATIC}"
  INTERFACE "${COLAMD_STATIC}"
  INTERFACE "${BTF_STATIC}"
  INTERFACE "${AMD_STATIC}"
  INTERFACE "${SUITESPARSE_CONFIG_STATIC}")
set_target_properties(SUNDIALS::KLU PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                               "${KLU_INCLUDE_DIR}")

find_dependency(SUNDIALS REQUIRED PATHS
             "@CMAKE_SOURCE_DIR@/ThirdParty/sundials/build/@CMAKE_INSTALL_LIBDIR@/cmake/sundials/")

if(@Boost_CHRONO_FOUND@)
  find_dependency(Boost COMPONENTS chrono REQUIRED)
endif()

if(@HDF5_FOUND@)
  find_dependency(
    HDF5
    COMPONENTS C HL CXX
    REQUIRED)
endif()

include("${CMAKE_CURRENT_LIST_DIR}/AmiciTargets.cmake")

check_required_components(Amici)
