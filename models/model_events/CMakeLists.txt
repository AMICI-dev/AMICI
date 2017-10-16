project(model_events)

cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wno-unused-function -Wno-unused-variable -Wno-unused-but-set-variable")

set(AMICI_DIR "${CMAKE_CURRENT_LIST_DIR}/../../")
set(MODEL_DIR "${AMICI_DIR}/models/model_events")
set(SUITESPARSE_DIR "${AMICI_DIR}/SuiteSparse/")
set(SUITESPARSE_LIB_DIR "${AMICI_DIR}/SuiteSparse/lib")
set(SUNDIALS_LIB_DIR "${AMICI_DIR}/sundials/build/lib")

find_package(HDF5 COMPONENTS C HL REQUIRED)

set(BLAS "CBLAS" CACHE STRING "BLAS library to use")
set_property(CACHE BLAS PROPERTY STRINGS "CBLAS" "MKL")

if(${BLAS} STREQUAL "MKL")
    set(BLAS_LIBRARIES -lmkl CACHE STRING "")
else()
    set(BLAS_LIBRARIES -lcblas CACHE STRING "")
endif()

include_directories("${AMICI_DIR}")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories("${HDF5_INCLUDE_DIRS}")
include_directories("${AMICI_DIR}/sundials/include")
include_directories("${SUITESPARSE_DIR}/include")

set(SRC_LIST_LIB ${MODEL_DIR}/model_events_xdot.cpp
${MODEL_DIR}/model_events_w.cpp
${MODEL_DIR}/model_events_dwdx.cpp
${MODEL_DIR}/model_events_J.cpp
${MODEL_DIR}/model_events_x0.cpp
${MODEL_DIR}/model_events_Jv.cpp
${MODEL_DIR}/model_events_JBand.cpp
${MODEL_DIR}/model_events_JSparse.cpp
${MODEL_DIR}/model_events_JDiag.cpp
${MODEL_DIR}/model_events_y.cpp
${MODEL_DIR}/model_events_z.cpp
${MODEL_DIR}/model_events_rz.cpp
${MODEL_DIR}/model_events_deltax.cpp
${MODEL_DIR}/model_events_root.cpp
${MODEL_DIR}/model_events_Jy.cpp
${MODEL_DIR}/model_events_Jz.cpp
${MODEL_DIR}/model_events_Jrz.cpp
${MODEL_DIR}/model_events_sigma_y.cpp
${MODEL_DIR}/model_events_sigma_z.cpp
${MODEL_DIR}/model_events_sxdot.cpp
${MODEL_DIR}/model_events_sx0.cpp
${MODEL_DIR}/model_events_sz.cpp
${MODEL_DIR}/model_events_deltasx.cpp
${MODEL_DIR}/model_events_stau.cpp
${MODEL_DIR}/model_events_srz.cpp
${MODEL_DIR}/model_events_dJydy.cpp
${MODEL_DIR}/model_events_dJydsigma.cpp
${MODEL_DIR}/model_events_dJzdz.cpp
${MODEL_DIR}/model_events_dJzdsigma.cpp
${MODEL_DIR}/model_events_dJrzdz.cpp
${MODEL_DIR}/model_events_dJrzdsigma.cpp
${MODEL_DIR}/model_events_dwdp.cpp
${MODEL_DIR}/model_events_dxdotdp.cpp
${MODEL_DIR}/model_events_dydp.cpp
${MODEL_DIR}/model_events_dsigma_ydp.cpp
${MODEL_DIR}/model_events_dsigma_zdp.cpp
${MODEL_DIR}/model_events_dydx.cpp
${MODEL_DIR}/model_events_dzdx.cpp
${MODEL_DIR}/model_events_dzdp.cpp
${MODEL_DIR}/model_events_drzdx.cpp
${MODEL_DIR}/model_events_drzdp.cpp
${MODEL_DIR}/model_events_xBdot.cpp
${MODEL_DIR}/model_events_qBdot.cpp
${MODEL_DIR}/model_events_JB.cpp
${MODEL_DIR}/model_events_JvB.cpp
${MODEL_DIR}/model_events_JBandB.cpp
${MODEL_DIR}/model_events_JSparseB.cpp
${MODEL_DIR}/model_events_deltaxB.cpp
${MODEL_DIR}/model_events_deltaqB.cpp

${MODEL_DIR}/wrapfunctions.cpp
)
    
add_library(${PROJECT_NAME} ${SRC_LIST_LIB})

if(APPLE)
    set(SHARED_OBJECT_EXTENSION .dylib)
else()
    set(SHARED_OBJECT_EXTENSION .so)
endif()
    
target_link_libraries(${PROJECT_NAME}
    ${AMICI_DIR}/build/libamici.a
    ${SUNDIALS_LIB_DIR}/libsundials_nvecserial${SHARED_OBJECT_EXTENSION}
    ${SUNDIALS_LIB_DIR}/libsundials_cvodes${SHARED_OBJECT_EXTENSION}
    ${SUNDIALS_LIB_DIR}/libsundials_idas${SHARED_OBJECT_EXTENSION}
    ${SUITESPARSE_LIB_DIR}/libcolamd${SHARED_OBJECT_EXTENSION}
    ${SUITESPARSE_LIB_DIR}/libklu${SHARED_OBJECT_EXTENSION}
    ${SUITESPARSE_LIB_DIR}/libbtf${SHARED_OBJECT_EXTENSION}
    ${SUITESPARSE_LIB_DIR}/libamd${SHARED_OBJECT_EXTENSION}
    ${SUITESPARSE_LIB_DIR}/libsuitesparseconfig${SHARED_OBJECT_EXTENSION}
    ${HDF5_HL_LIBRARIES}
    ${HDF5_C_LIBRARIES}
    ${BLAS_LIBRARIES}
    -ldl -lz -lm
)

set(SRC_LIST_EXE main.cpp)

add_executable(simulate_${PROJECT_NAME} ${SRC_LIST_EXE})
    
target_link_libraries(simulate_${PROJECT_NAME} ${AMICI_DIR}/build/libamici.a ${PROJECT_NAME})

if($ENV{ENABLE_GCOV_COVERAGE})
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 --coverage")
endif()
