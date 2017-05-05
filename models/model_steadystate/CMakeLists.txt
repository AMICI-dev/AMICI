project(model_steadystate)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wno-unused-function")
add_definitions(-DAMICI_WITHOUT_MATLAB)

set(AMICI_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../")
set(MODEL_DIR "${AMICI_DIR}/models/model_steadystate")
set(SUITESPARSE_DIR "${AMICI_DIR}/SuiteSparse/")
set(SUITESPARSE_LIB_DIR "${AMICI_DIR}/SuiteSparse/lib")
set(SUNDIALS_LIB_DIR "${AMICI_DIR}/sundials/build/lib")

find_package(HDF5 COMPONENTS C HL REQUIRED)
include_directories("${AMICI_DIR}")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories("${MODEL_DIR}")
include_directories("${HDF5_INCLUDE_DIRS}")
include_directories("${AMICI_DIR}/sundials/include")
include_directories("${SUITESPARSE_DIR}/include")


set(SRC_LIST_LIB
${MODEL_DIR}/wrapfunctions.cpp
${AMICI_DIR}/src/symbolic_functions.cpp
${AMICI_DIR}/src/amici.cpp
${AMICI_DIR}/src/udata.cpp
${AMICI_DIR}/src/rdata.cpp
${AMICI_DIR}/src/edata.cpp
${AMICI_DIR}/src/ami_hdf5.cpp
${AMICI_DIR}/src/spline.cpp
${MODEL_DIR}/model_steadystate_J.cpp
${MODEL_DIR}/model_steadystate_JB.cpp
${MODEL_DIR}/model_steadystate_JBand.cpp
${MODEL_DIR}/model_steadystate_JBandB.cpp
${MODEL_DIR}/model_steadystate_JSparse.cpp
${MODEL_DIR}/model_steadystate_JSparseB.cpp
${MODEL_DIR}/model_steadystate_Jv.cpp
${MODEL_DIR}/model_steadystate_JvB.cpp
${MODEL_DIR}/model_steadystate_Jy.cpp
${MODEL_DIR}/model_steadystate_Jz.cpp
${MODEL_DIR}/model_steadystate_dJydp.cpp
${MODEL_DIR}/model_steadystate_dJydx.cpp
${MODEL_DIR}/model_steadystate_dJydy.cpp
${MODEL_DIR}/model_steadystate_dJzdp.cpp
${MODEL_DIR}/model_steadystate_dJzdx.cpp
${MODEL_DIR}/model_steadystate_deltaqB.cpp
${MODEL_DIR}/model_steadystate_deltasx.cpp
${MODEL_DIR}/model_steadystate_deltax.cpp
${MODEL_DIR}/model_steadystate_deltaxB.cpp
${MODEL_DIR}/model_steadystate_dsigma_ydp.cpp
${MODEL_DIR}/model_steadystate_dsigma_zdp.cpp
${MODEL_DIR}/model_steadystate_dwdp.cpp
${MODEL_DIR}/model_steadystate_dwdx.cpp
${MODEL_DIR}/model_steadystate_dxdotdp.cpp
${MODEL_DIR}/model_steadystate_dydp.cpp
${MODEL_DIR}/model_steadystate_dydx.cpp
${MODEL_DIR}/model_steadystate_dzdp.cpp
${MODEL_DIR}/model_steadystate_dzdx.cpp
${MODEL_DIR}/model_steadystate_qBdot.cpp
${MODEL_DIR}/model_steadystate_root.cpp
${MODEL_DIR}/model_steadystate_sJy.cpp
${MODEL_DIR}/model_steadystate_sJz.cpp
${MODEL_DIR}/model_steadystate_sigma_y.cpp
${MODEL_DIR}/model_steadystate_sigma_z.cpp
${MODEL_DIR}/model_steadystate_sroot.cpp
${MODEL_DIR}/model_steadystate_stau.cpp
${MODEL_DIR}/model_steadystate_sx0.cpp
${MODEL_DIR}/model_steadystate_sxdot.cpp
${MODEL_DIR}/model_steadystate_sy.cpp
${MODEL_DIR}/model_steadystate_sz.cpp
${MODEL_DIR}/model_steadystate_sz_tf.cpp
${MODEL_DIR}/model_steadystate_w.cpp
${MODEL_DIR}/model_steadystate_x0.cpp
${MODEL_DIR}/model_steadystate_xBdot.cpp
${MODEL_DIR}/model_steadystate_xdot.cpp
${MODEL_DIR}/model_steadystate_y.cpp
${MODEL_DIR}/model_steadystate_z.cpp
)

add_library(${PROJECT_NAME} ${SRC_LIST_LIB})

if(APPLE)
    set(SHARED_OBJECT_EXTENSION .dylib)
else()
    set(SHARED_OBJECT_EXTENSION .so)
endif()

target_link_libraries(${PROJECT_NAME}
"${SUNDIALS_LIB_DIR}/libsundials_nvecserial${SHARED_OBJECT_EXTENSION}"
"${SUNDIALS_LIB_DIR}/libsundials_cvodes${SHARED_OBJECT_EXTENSION}"
"${SUITESPARSE_LIB_DIR}/libcolamd${SHARED_OBJECT_EXTENSION}"
"${SUITESPARSE_LIB_DIR}/libklu${SHARED_OBJECT_EXTENSION}"
"${SUITESPARSE_LIB_DIR}/libbtf${SHARED_OBJECT_EXTENSION}"
"${SUITESPARSE_LIB_DIR}/libamd${SHARED_OBJECT_EXTENSION}"
"${SUITESPARSE_LIB_DIR}/libsuitesparseconfig${SHARED_OBJECT_EXTENSION}"
"${HDF5_HL_LIBRARIES}"
"${HDF5_C_LIBRARIES}"
"-lpthread -ldl -lz"
"-lm"
)

set(SRC_LIST_EXE main.cpp)

add_executable(simulate_${PROJECT_NAME} ${SRC_LIST_EXE})

target_link_libraries(simulate_${PROJECT_NAME} ${PROJECT_NAME})

