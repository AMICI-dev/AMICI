@PACKAGE_INIT@

include(CMakeFindDependencyMacro)

find_package(OpenMP)

# Current SUNDIALSConfig.cmake doesn't use KLU's FindKLU, but hardcodes paths,
# and is not relocatable. This does not work with Python package installation in
# tmpdirs.
find_package(SuiteSparse_config REQUIRED)
find_package(AMD REQUIRED)
find_package(BTF REQUIRED)
find_package(COLAMD REQUIRED)
find_package(KLU REQUIRED)

add_library(SUNDIALS::KLU INTERFACE IMPORTED)
target_link_libraries(SUNDIALS::KLU
    INTERFACE "${KLU_STATIC}"
    INTERFACE "${COLAMD_STATIC}"
    INTERFACE "${BTF_STATIC}"
    INTERFACE "${AMD_STATIC}"
    INTERFACE "${SUITESPARSE_CONFIG_STATIC}"
)
set_target_properties(
    SUNDIALS::KLU PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${KLU_INCLUDE_DIR}"
)


find_package(SUNDIALS REQUIRED PATHS
             "@CMAKE_SOURCE_DIR@/ThirdParty/sundials/build/lib/cmake/sundials/")

if(@HDF5_FOUND@)
    find_package(
      HDF5
      COMPONENTS C HL CXX
      REQUIRED)
endif()

include("${CMAKE_CURRENT_LIST_DIR}/AmiciTargets.cmake")

check_required_components(Amici)
