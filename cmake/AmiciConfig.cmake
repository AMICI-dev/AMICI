include("${CMAKE_CURRENT_LIST_DIR}/AmiciTargets.cmake")
set_property(TARGET Upstream::amici PROPERTY INCLUDE_DIRECTORIES @AmiciConfigIncludes@)
