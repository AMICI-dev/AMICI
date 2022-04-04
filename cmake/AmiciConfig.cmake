@PACKAGE_INIT@

include(CMakeFindDependencyMacro)

find_package(SUNDIALS REQUIRED PATHS "@CMAKE_SOURCE_DIR@/ThirdParty/sundials/build/lib/cmake/sundials/")

include("${CMAKE_CURRENT_LIST_DIR}/AmiciTargets.cmake")

check_required_components(Amici)
