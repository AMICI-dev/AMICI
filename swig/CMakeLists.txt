#
# Build additional AMICI interfaces via swig
#

find_package(SWIG REQUIRED 3)
include(${SWIG_USE_FILE})

set(AMICI_INTERFACE_LIST
    ${CMAKE_CURRENT_SOURCE_DIR}/amici.i
    ${CMAKE_CURRENT_SOURCE_DIR}/edata.i
    ${CMAKE_CURRENT_SOURCE_DIR}/rdata.i
    ${CMAKE_CURRENT_SOURCE_DIR}/model.i
    ${CMAKE_CURRENT_SOURCE_DIR}/model_ode.i
    ${CMAKE_CURRENT_SOURCE_DIR}/model_dae.i
    ${CMAKE_CURRENT_SOURCE_DIR}/solver.i
    ${CMAKE_CURRENT_SOURCE_DIR}/solver_cvodes.i
    ${CMAKE_CURRENT_SOURCE_DIR}/solver_idas.i
    ${CMAKE_CURRENT_SOURCE_DIR}/std_unique_ptr.i
    ${CMAKE_CURRENT_SOURCE_DIR}/hdf5.i
    ${CMAKE_CURRENT_SOURCE_DIR}/abstract_model.i
    ${CMAKE_CURRENT_SOURCE_DIR}/stdvec2numpy.h
)

# Add target to show files in IDE
add_custom_target(swigInterface SOURCES ${AMICI_INTERFACE_LIST})
