project(model_dirac)
cmake_minimum_required(VERSION 2.8)
set(cmake_build_type Debug)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wno-unused-function")
add_definitions(-DAMICI_WITHOUT_MATLAB)

set(AMICI_DIR "${CMAKE_SOURCE_DIR}")
set(MODEL_DIR "${AMICI_DIR}/models/model_dirac")
set(SUITESPARSE_DIR "${AMICI_DIR}/SuiteSparse/")
set(SUITESPARSE_LIB_DIR "${AMICI_DIR}/SuiteSparse/lib")
set(SUNDIALS_LIB_DIR "${AMICI_DIR}/build/sundials/lib")

set(SRC_LIST
${MODEL_DIR}/main.cpp
${MODEL_DIR}/wrapfunctions.cpp
${AMICI_DIR}/src/symbolic_functions.cpp
${AMICI_DIR}/src/amici.cpp
${AMICI_DIR}/src/udata.cpp
${AMICI_DIR}/src/rdata.cpp
${AMICI_DIR}/src/edata.cpp
${AMICI_DIR}/src/ami_hdf5.cpp
${MODEL_DIR}/model_dirac_J.cpp
${MODEL_DIR}/model_dirac_JB.cpp
${MODEL_DIR}/model_dirac_JBand.cpp
${MODEL_DIR}/model_dirac_JBandB.cpp
${MODEL_DIR}/model_dirac_JSparse.cpp
${MODEL_DIR}/model_dirac_JSparseB.cpp
${MODEL_DIR}/model_dirac_Jv.cpp
${MODEL_DIR}/model_dirac_JvB.cpp
${MODEL_DIR}/model_dirac_Jy.cpp
${MODEL_DIR}/model_dirac_Jz.cpp
${MODEL_DIR}/model_dirac_dJydp.cpp
${MODEL_DIR}/model_dirac_dJydx.cpp
${MODEL_DIR}/model_dirac_dJzdp.cpp
${MODEL_DIR}/model_dirac_dJzdx.cpp
${MODEL_DIR}/model_dirac_deltaqB.cpp
${MODEL_DIR}/model_dirac_deltasx.cpp
${MODEL_DIR}/model_dirac_deltax.cpp
${MODEL_DIR}/model_dirac_deltaxB.cpp
${MODEL_DIR}/model_dirac_dsigma_ydp.cpp
${MODEL_DIR}/model_dirac_dsigma_zdp.cpp
${MODEL_DIR}/model_dirac_dwdp.cpp
${MODEL_DIR}/model_dirac_dwdx.cpp
${MODEL_DIR}/model_dirac_dxdotdp.cpp
${MODEL_DIR}/model_dirac_dydp.cpp
${MODEL_DIR}/model_dirac_dydx.cpp
${MODEL_DIR}/model_dirac_dzdp.cpp
${MODEL_DIR}/model_dirac_dzdx.cpp
${MODEL_DIR}/model_dirac_qBdot.cpp
${MODEL_DIR}/model_dirac_root.cpp
${MODEL_DIR}/model_dirac_sJy.cpp
${MODEL_DIR}/model_dirac_sJz.cpp
${MODEL_DIR}/model_dirac_sigma_y.cpp
${MODEL_DIR}/model_dirac_sigma_z.cpp
${MODEL_DIR}/model_dirac_sroot.cpp
${MODEL_DIR}/model_dirac_stau.cpp
${MODEL_DIR}/model_dirac_sx0.cpp
${MODEL_DIR}/model_dirac_sxdot.cpp
${MODEL_DIR}/model_dirac_sy.cpp
${MODEL_DIR}/model_dirac_sz.cpp
${MODEL_DIR}/model_dirac_sz_tf.cpp
${MODEL_DIR}/model_dirac_w.cpp
${MODEL_DIR}/model_dirac_x0.cpp
${MODEL_DIR}/model_dirac_xBdot.cpp
${MODEL_DIR}/model_dirac_xdot.cpp
${MODEL_DIR}/model_dirac_y.cpp
${MODEL_DIR}/model_dirac_z.cpp
)

include_directories("${AMICI_DIR}")
include_directories("${AMICI_DIR}/models/model_dirac")
include_directories("${AMICI_DIR}/sundials/include")
include_directories("${SUITESPARSE_DIR}/KLU/Include")
include_directories("${SUITESPARSE_DIR}/AMD/Include")
include_directories("${SUITESPARSE_DIR}/SuiteSparse_config")
include_directories("${SUITESPARSE_DIR}/COLAMD/Include")
include_directories("${SUITESPARSE_DIR}/BTF/Include")

add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}
"-lpthread -ldl -lz -lhdf5_hl -lhdf5 -L${SUNDIALS_LIB_DIR} -lsundials_cvodes -lsundials_nvecserial -L${SUITESPARSE_LIB_DIR}  -lklu -lbtf -lamd -lcolamd -lsuitesparseconfig"
"-lm"
)
