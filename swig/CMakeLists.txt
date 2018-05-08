#
# Build additional AMICI interfaces via swig
#

#find_package(SWIG REQUIRED 3)
#include(${SWIG_USE_FILE})

#set(AMICI_INTERFACE_LIST
#    ${AMICI_DIR}/swig/edata.i
#    ${AMICI_DIR}/swig/rdata.i
#    ${AMICI_DIR}/swig/model.i
#    ${AMICI_DIR}/swig/model_ode.i
#    ${AMICI_DIR}/swig/model_dae.i
#    ${AMICI_DIR}/swig/solver.i
#    ${AMICI_DIR}/swig/solver_cvodes.i
#    ${AMICI_DIR}/swig/solver_idas.i
#    ${AMICI_DIR}/swig/std_unique_ptr.i
#)

## Add subdirectories for each language if desired
#option(BUILD_PYTHON "Build Python SWIG module" ON)
#if(BUILD_PYTHON)
#    add_subdirectory(python)
#endif()

