project(model_jakstat_adjoint)

cmake_minimum_required(VERSION 2.8)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wno-unused-function -Wno-unused-variable")
add_definitions(-DAMICI_WITHOUT_MATLAB)
    
set(AMICI_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../")
set(MODEL_DIR "${AMICI_DIR}/models/model_jakstat_adjoint")
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
include_directories("${MODEL_DIR}")
include_directories("${HDF5_INCLUDE_DIRS}")
include_directories("${AMICI_DIR}/sundials/include")
include_directories("${SUITESPARSE_DIR}/include")


set(SRC_LIST_LIB ${MODEL_DIR}/model_jakstat_adjoint_J.cpp
${MODEL_DIR}/model_jakstat_adjoint_JB.cpp
${MODEL_DIR}/model_jakstat_adjoint_JBand.cpp
${MODEL_DIR}/model_jakstat_adjoint_JBandB.cpp
${MODEL_DIR}/model_jakstat_adjoint_JSparse.cpp
${MODEL_DIR}/model_jakstat_adjoint_JSparseB.cpp
${MODEL_DIR}/model_jakstat_adjoint_Jv.cpp
${MODEL_DIR}/model_jakstat_adjoint_JvB.cpp
${MODEL_DIR}/model_jakstat_adjoint_Jy.cpp
${MODEL_DIR}/model_jakstat_adjoint_Jz.cpp
${MODEL_DIR}/model_jakstat_adjoint_dJydsigma.cpp
${MODEL_DIR}/model_jakstat_adjoint_dJydy.cpp
${MODEL_DIR}/model_jakstat_adjoint_dJzdsigma.cpp
${MODEL_DIR}/model_jakstat_adjoint_dJzdz.cpp
${MODEL_DIR}/model_jakstat_adjoint_deltaqB.cpp
${MODEL_DIR}/model_jakstat_adjoint_deltasx.cpp
${MODEL_DIR}/model_jakstat_adjoint_deltax.cpp
${MODEL_DIR}/model_jakstat_adjoint_deltaxB.cpp
${MODEL_DIR}/model_jakstat_adjoint_drootdp.cpp
${MODEL_DIR}/model_jakstat_adjoint_drootdx.cpp
${MODEL_DIR}/model_jakstat_adjoint_dsigma_ydp.cpp
${MODEL_DIR}/model_jakstat_adjoint_dsigma_zdp.cpp
${MODEL_DIR}/model_jakstat_adjoint_dwdp.cpp
${MODEL_DIR}/model_jakstat_adjoint_dwdx.cpp
${MODEL_DIR}/model_jakstat_adjoint_dxdotdp.cpp
${MODEL_DIR}/model_jakstat_adjoint_dydp.cpp
${MODEL_DIR}/model_jakstat_adjoint_dydx.cpp
${MODEL_DIR}/model_jakstat_adjoint_dzdp.cpp
${MODEL_DIR}/model_jakstat_adjoint_dzdx.cpp
${MODEL_DIR}/model_jakstat_adjoint_qBdot.cpp
${MODEL_DIR}/model_jakstat_adjoint_root.cpp
${MODEL_DIR}/model_jakstat_adjoint_sigma_y.cpp
${MODEL_DIR}/model_jakstat_adjoint_sigma_z.cpp
${MODEL_DIR}/model_jakstat_adjoint_sroot.cpp
${MODEL_DIR}/model_jakstat_adjoint_stau.cpp
${MODEL_DIR}/model_jakstat_adjoint_sx0.cpp
${MODEL_DIR}/model_jakstat_adjoint_sxdot.cpp
${MODEL_DIR}/model_jakstat_adjoint_sz.cpp
${MODEL_DIR}/model_jakstat_adjoint_w.cpp
${MODEL_DIR}/model_jakstat_adjoint_x0.cpp
${MODEL_DIR}/model_jakstat_adjoint_xBdot.cpp
${MODEL_DIR}/model_jakstat_adjoint_xdot.cpp
${MODEL_DIR}/model_jakstat_adjoint_y.cpp
${MODEL_DIR}/model_jakstat_adjoint_z.cpp
${MODEL_DIR}/wrapfunctions.cpp
${AMICI_DIR}/src/symbolic_functions.cpp
${AMICI_DIR}/src/amici_interface_cpp.cpp
${AMICI_DIR}/src/amici.cpp
${AMICI_DIR}/src/udata.cpp
${AMICI_DIR}/src/rdata.cpp
${AMICI_DIR}/src/tdata.cpp
${AMICI_DIR}/src/edata.cpp
${AMICI_DIR}/src/ami_hdf5.cpp
${AMICI_DIR}/src/spline.cpp)
    
add_library(${PROJECT_NAME} ${SRC_LIST_LIB})

if(APPLE)
    set(SHARED_OBJECT_EXTENSION .dylib)
else()
    set(SHARED_OBJECT_EXTENSION .so)
endif()
    
target_link_libraries(${PROJECT_NAME} ${SUNDIALS_LIB_DIR}/libsundials_nvecserial${SHARED_OBJECT_EXTENSION}
${SUNDIALS_LIB_DIR}/libsundials_cvodes${SHARED_OBJECT_EXTENSION}
${SUITESPARSE_LIB_DIR}/libcolamd${SHARED_OBJECT_EXTENSION}
${SUITESPARSE_LIB_DIR}/libklu${SHARED_OBJECT_EXTENSION}
${SUITESPARSE_LIB_DIR}/libbtf${SHARED_OBJECT_EXTENSION}
${SUITESPARSE_LIB_DIR}/libamd${SHARED_OBJECT_EXTENSION}
${SUITESPARSE_LIB_DIR}/libsuitesparseconfig${SHARED_OBJECT_EXTENSION}
${HDF5_HL_LIBRARIES}
${HDF5_C_LIBRARIES}
-lpthread -ldl -lz -lcblas
-lm)

set(SRC_LIST_EXE main.cpp)

add_executable(simulate_${PROJECT_NAME} ${SRC_LIST_EXE})
    
target_link_libraries(simulate_${PROJECT_NAME} ${PROJECT_NAME})
