#
# Build additional AMICI interfaces via swig
#

find_package(SWIG REQUIRED 3)
include(${SWIG_USE_FILE})

set(AMICI_INTERFACE_LIST
    ${CMAKE_CURRENT_LIST_DIR}/amici.i
    ${CMAKE_CURRENT_LIST_DIR}/edata.i
    ${CMAKE_CURRENT_LIST_DIR}/rdata.i
    ${CMAKE_CURRENT_LIST_DIR}/model.i
    ${CMAKE_CURRENT_LIST_DIR}/model_ode.i
    ${CMAKE_CURRENT_LIST_DIR}/model_dae.i
    ${CMAKE_CURRENT_LIST_DIR}/solver.i
    ${CMAKE_CURRENT_LIST_DIR}/solver_cvodes.i
    ${CMAKE_CURRENT_LIST_DIR}/solver_idas.i
    ${CMAKE_CURRENT_LIST_DIR}/std_unique_ptr.i
)

# Add target to show files in IDE
add_custom_target(swigInterface SOURCES ${AMICI_INTERFACE_LIST})

## Add subdirectories for each language if desired
#option(BUILD_PYTHON "Build Python SWIG module" ON)
#if(BUILD_PYTHON)
#    add_subdirectory(python)
#endif()

