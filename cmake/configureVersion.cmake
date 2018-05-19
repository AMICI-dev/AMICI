include(${CMAKE_CURRENT_LIST_DIR}/version.cmake)
set(AMICI_VERSION ${PROJECT_VERSION})
configure_file(${SRC} ${DST} @ONLY)


