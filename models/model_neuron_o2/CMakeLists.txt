project(model_neuron_o2)

cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wno-unused-function -Wno-unused-variable -Wno-unused-but-set-variable")

set(AMICI_DIR "${CMAKE_CURRENT_LIST_DIR}/../../")
set(MODEL_DIR "${AMICI_DIR}/models/model_neuron_o2")
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

set(SRC_LIST_LIB ${MODEL_DIR}/model_neuron_o2_xdot.cpp
${MODEL_DIR}/model_neuron_o2_w.cpp
${MODEL_DIR}/model_neuron_o2_dwdx.cpp
${MODEL_DIR}/model_neuron_o2_J.cpp
${MODEL_DIR}/model_neuron_o2_x0.cpp
${MODEL_DIR}/model_neuron_o2_Jv.cpp
${MODEL_DIR}/model_neuron_o2_JBand.cpp
${MODEL_DIR}/model_neuron_o2_JSparse.cpp
${MODEL_DIR}/model_neuron_o2_JDiag.cpp
${MODEL_DIR}/model_neuron_o2_y.cpp
${MODEL_DIR}/model_neuron_o2_z.cpp
${MODEL_DIR}/model_neuron_o2_rz.cpp
${MODEL_DIR}/model_neuron_o2_deltax.cpp
${MODEL_DIR}/model_neuron_o2_root.cpp
${MODEL_DIR}/model_neuron_o2_Jy.cpp
${MODEL_DIR}/model_neuron_o2_Jz.cpp
${MODEL_DIR}/model_neuron_o2_Jrz.cpp
${MODEL_DIR}/model_neuron_o2_sigma_y.cpp
${MODEL_DIR}/model_neuron_o2_sigma_z.cpp
${MODEL_DIR}/model_neuron_o2_sxdot.cpp
${MODEL_DIR}/model_neuron_o2_sx0.cpp
${MODEL_DIR}/model_neuron_o2_sz.cpp
${MODEL_DIR}/model_neuron_o2_deltasx.cpp
${MODEL_DIR}/model_neuron_o2_stau.cpp
${MODEL_DIR}/model_neuron_o2_srz.cpp
${MODEL_DIR}/model_neuron_o2_dJydy.cpp
${MODEL_DIR}/model_neuron_o2_dJydsigma.cpp
${MODEL_DIR}/model_neuron_o2_dJzdz.cpp
${MODEL_DIR}/model_neuron_o2_dJzdsigma.cpp
${MODEL_DIR}/model_neuron_o2_dJrzdz.cpp
${MODEL_DIR}/model_neuron_o2_dJrzdsigma.cpp
${MODEL_DIR}/model_neuron_o2_dwdp.cpp
${MODEL_DIR}/model_neuron_o2_dxdotdp.cpp
${MODEL_DIR}/model_neuron_o2_dydp.cpp
${MODEL_DIR}/model_neuron_o2_dsigma_ydp.cpp
${MODEL_DIR}/model_neuron_o2_dsigma_zdp.cpp
${MODEL_DIR}/model_neuron_o2_dydx.cpp
${MODEL_DIR}/model_neuron_o2_dzdx.cpp
${MODEL_DIR}/model_neuron_o2_dzdp.cpp
${MODEL_DIR}/model_neuron_o2_drzdx.cpp
${MODEL_DIR}/model_neuron_o2_drzdp.cpp

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
