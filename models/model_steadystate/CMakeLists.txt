project(model_steadystate)

cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wno-unused-function -Wno-unused-variable -Wno-unused-but-set-variable")

add_definitions(-DAMICI_WITHOUT_MATLAB)
    
set(AMICI_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../")
set(MODEL_DIR "${AMICI_DIR}/models/model_steadystate")
set(SUITESPARSE_DIR "${AMICI_DIR}/SuiteSparse/")
set(SUITESPARSE_LIB_DIR "${AMICI_DIR}/SuiteSparse/lib")
set(SUNDIALS_LIB_DIR "${AMICI_DIR}/sundials/build/lib")

find_package(HDF5 COMPONENTS C HL REQUIRED)

set(BLAS "CBLAS" CACHE STRING "BLAS library to use")
set_property(CACHE BLAS PROPERTY STRINGS "CBLAS" "MKL")

if(${BLAS} STREQUAL "MKL")
    add_definitions(-DAMICI_BLAS_MKL)
    set(BLAS_LIBRARIES -lmkl CACHE STRING "")
    set(BLAS_INCLUDE_DIRS "" CACHE STRING "")
else()
    add_definitions(-DAMICI_BLAS_CBLAS)
    set(BLAS_INCLUDE_DIRS "" CACHE STRING "")
    set(BLAS_LIBRARIES -lcblas CACHE STRING "")
endif()

include_directories("${AMICI_DIR}")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories("${HDF5_INCLUDE_DIRS}")
include_directories("${AMICI_DIR}/sundials/include")
include_directories("${SUITESPARSE_DIR}/include")

if(NOT "${BLAS_INCLUDE_DIRS}" STREQUAL "")
    include_directories(${BLAS_INCLUDE_DIRS})
endif()

set(SRC_LIST_LIB ${MODEL_DIR}/model_steadystate_J.cpp
${MODEL_DIR}/model_steadystate_JB.cpp
${MODEL_DIR}/model_steadystate_JBand.cpp
${MODEL_DIR}/model_steadystate_JBandB.cpp
${MODEL_DIR}/model_steadystate_JSparse.cpp
${MODEL_DIR}/model_steadystate_JSparseB.cpp
${MODEL_DIR}/model_steadystate_Jrz.cpp
${MODEL_DIR}/model_steadystate_Jv.cpp
${MODEL_DIR}/model_steadystate_JvB.cpp
${MODEL_DIR}/model_steadystate_Jy.cpp
${MODEL_DIR}/model_steadystate_Jz.cpp
${MODEL_DIR}/model_steadystate_dJrzdsigma.cpp
${MODEL_DIR}/model_steadystate_dJrzdz.cpp
${MODEL_DIR}/model_steadystate_dJydsigma.cpp
${MODEL_DIR}/model_steadystate_dJydy.cpp
${MODEL_DIR}/model_steadystate_dJzdsigma.cpp
${MODEL_DIR}/model_steadystate_dJzdz.cpp
${MODEL_DIR}/model_steadystate_deltaqB.cpp
${MODEL_DIR}/model_steadystate_deltasx.cpp
${MODEL_DIR}/model_steadystate_deltax.cpp
${MODEL_DIR}/model_steadystate_deltaxB.cpp
${MODEL_DIR}/model_steadystate_drzdp.cpp
${MODEL_DIR}/model_steadystate_drzdx.cpp
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
${MODEL_DIR}/model_steadystate_rz.cpp
${MODEL_DIR}/model_steadystate_sigma_y.cpp
${MODEL_DIR}/model_steadystate_sigma_z.cpp
${MODEL_DIR}/model_steadystate_srz.cpp
${MODEL_DIR}/model_steadystate_stau.cpp
${MODEL_DIR}/model_steadystate_sx0.cpp
${MODEL_DIR}/model_steadystate_sxdot.cpp
${MODEL_DIR}/model_steadystate_sz.cpp
${MODEL_DIR}/model_steadystate_w.cpp
${MODEL_DIR}/model_steadystate_x0.cpp
${MODEL_DIR}/model_steadystate_xBdot.cpp
${MODEL_DIR}/model_steadystate_xdot.cpp
${MODEL_DIR}/model_steadystate_y.cpp
${MODEL_DIR}/model_steadystate_z.cpp

${MODEL_DIR}/wrapfunctions.cpp
${AMICI_DIR}/src/symbolic_functions.cpp
${AMICI_DIR}/src/amici_interface_cpp.cpp
${AMICI_DIR}/src/amici.cpp
${AMICI_DIR}/src/amici_misc.cpp
${AMICI_DIR}/src/udata.cpp
${AMICI_DIR}/src/rdata.cpp
${AMICI_DIR}/src/tdata.cpp
${AMICI_DIR}/src/edata.cpp
${AMICI_DIR}/src/amici_hdf5.cpp
${AMICI_DIR}/src/spline.cpp)
    
add_library(${PROJECT_NAME} ${SRC_LIST_LIB})

if(APPLE)
    set(SHARED_OBJECT_EXTENSION .dylib)
else()
    set(SHARED_OBJECT_EXTENSION .so)
endif()
    
target_link_libraries(${PROJECT_NAME}
    ${SUNDIALS_LIB_DIR}/libsundials_nvecserial${SHARED_OBJECT_EXTENSION}
    ${SUNDIALS_LIB_DIR}/libsundials_cvodes${SHARED_OBJECT_EXTENSION}
    ${SUITESPARSE_LIB_DIR}/libcolamd${SHARED_OBJECT_EXTENSION}
    ${SUITESPARSE_LIB_DIR}/libklu${SHARED_OBJECT_EXTENSION}
    ${SUITESPARSE_LIB_DIR}/libbtf${SHARED_OBJECT_EXTENSION}
    ${SUITESPARSE_LIB_DIR}/libamd${SHARED_OBJECT_EXTENSION}
    ${SUITESPARSE_LIB_DIR}/libsuitesparseconfig${SHARED_OBJECT_EXTENSION}
    ${HDF5_HL_LIBRARIES}
    ${HDF5_C_LIBRARIES}
    ${BLAS_LIBRARIES}
    -lpthread -ldl -lz -lm
)

set(SRC_LIST_EXE main.cpp)

add_executable(simulate_${PROJECT_NAME} ${SRC_LIST_EXE})
    
target_link_libraries(simulate_${PROJECT_NAME} ${PROJECT_NAME})
