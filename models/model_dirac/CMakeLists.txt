project(model_dirac)
cmake_minimum_required(VERSION 2.8)

set(cmake_build_type Debug)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wno-unused-function")
add_definitions(-DAMICI_WITHOUT_MATLAB)

set(AMICI_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../")
set(SUITESPARSE_DIR "${AMICI_DIR}/SuiteSparse/")
set(SUITESPARSE_LIB_DIR "${AMICI_DIR}/SuiteSparse/")
set(SUNDIALS_LIB_DIR "${AMICI_DIR}/build/sundials/lib")

find_package(HDF5 COMPONENTS C HL REQUIRED)

set(SRC_LIST
main.cpp
wrapfunctions.cpp
${AMICI_DIR}/src/symbolic_functions.cpp
${AMICI_DIR}/src/amici.cpp
${AMICI_DIR}/src/udata.cpp
${AMICI_DIR}/src/rdata.cpp
${AMICI_DIR}/src/edata.cpp
${AMICI_DIR}/src/ami_hdf5.cpp
${AMICI_DIR}/src/spline.cpp
model_dirac_J.cpp
model_dirac_JB.cpp
model_dirac_JBand.cpp
model_dirac_JBandB.cpp
model_dirac_JSparse.cpp
model_dirac_JSparseB.cpp
model_dirac_Jv.cpp
model_dirac_JvB.cpp
model_dirac_Jy.cpp
model_dirac_Jz.cpp
model_dirac_dJydp.cpp
model_dirac_dJydx.cpp
model_dirac_dJydy.cpp
model_dirac_dJzdp.cpp
model_dirac_dJzdx.cpp
model_dirac_deltaqB.cpp
model_dirac_deltasx.cpp
model_dirac_deltax.cpp
model_dirac_deltaxB.cpp
model_dirac_dsigma_ydp.cpp
model_dirac_dsigma_zdp.cpp
model_dirac_dwdp.cpp
model_dirac_dwdx.cpp
model_dirac_dxdotdp.cpp
model_dirac_dydp.cpp
model_dirac_dydx.cpp
model_dirac_dzdp.cpp
model_dirac_dzdx.cpp
model_dirac_qBdot.cpp
model_dirac_root.cpp
model_dirac_sJy.cpp
model_dirac_sJz.cpp
model_dirac_sigma_y.cpp
model_dirac_sigma_z.cpp
model_dirac_sroot.cpp
model_dirac_stau.cpp
model_dirac_sx0.cpp
model_dirac_sxdot.cpp
model_dirac_sy.cpp
model_dirac_sz.cpp
model_dirac_sz_tf.cpp
model_dirac_w.cpp
model_dirac_x0.cpp
model_dirac_xBdot.cpp
model_dirac_xdot.cpp
model_dirac_y.cpp
model_dirac_z.cpp
)

include_directories("${AMICI_DIR}")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories("${HDF5_INCLUDE_DIRS}")
include_directories("${AMICI_DIR}/sundials/include")
include_directories("${SUITESPARSE_DIR}/KLU/Include")
include_directories("${SUITESPARSE_DIR}/AMD/Include")
include_directories("${SUITESPARSE_DIR}/SuiteSparse_config")
include_directories("${SUITESPARSE_DIR}/COLAMD/Include")
include_directories("${SUITESPARSE_DIR}/BTF/Include")

add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}
"${SUNDIALS_LIB_DIR}/libsundials_nvecserial.so"
"${SUNDIALS_LIB_DIR}/libsundials_cvodes.so"
"${SUITESPARSE_LIB_DIR}/lib/libcolamd.so"
"${SUITESPARSE_LIB_DIR}/KLU/Lib/libklu.a"
"${SUITESPARSE_LIB_DIR}/BTF/Lib/libbtf.a"
"${SUITESPARSE_LIB_DIR}/AMD/Lib/libamd.a"
"${SUITESPARSE_LIB_DIR}/COLAMD/Lib/libcolamd.a"
"${SUITESPARSE_LIB_DIR}/SuiteSparse_config/libsuitesparseconfig.a"
"${HDF5_HL_LIBRARIES}"
"${HDF5_C_LIBRARIES}"
"-lpthread -ldl -lz -lsz"
"-lm"
)
