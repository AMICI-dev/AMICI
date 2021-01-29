@PACKAGE_INIT@

include(CMakeFindDependencyMacro)

include("${CMAKE_CURRENT_LIST_DIR}/AmiciTargets.cmake")

check_required_components(Amici)
