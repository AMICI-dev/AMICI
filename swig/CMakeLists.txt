#
# Build additional AMICI interfaces via swig
#

# Use most recent SWIG version available
SET(CMAKE_FIND_PACKAGE_SORT_DIRECTION DEC)
find_package(SWIG REQUIRED)

set(SWIG_VERSION_MIN "3.0")
if (${SWIG_VERSION} VERSION_LESS ${SWIG_VERSION_MIN})
    message(FATAL_ERROR "Requiring SWIG>=${SWIG_VERSION_MIN} "
                        "but found only ${SWIG_VERSION}.")
endif()

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
