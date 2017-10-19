project(model_dirac)

cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wno-unused-function -Wno-unused-variable -Wno-unused-but-set-variable")

find_package(Amici HINTS ${CMAKE_CURRENT_LIST_DIR}/../../build)
get_property(AmiciConfigIncludes TARGET Upstream::amici PROPERTY INCLUDE_DIRECTORIES)

set(MODEL_DIR "${Amici_DIR}/../models/model_dirac")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories(${AmiciConfigIncludes})

set(SRC_LIST_LIB ${MODEL_DIR}/model_dirac_J.cpp
${MODEL_DIR}/model_dirac_JB.cpp
${MODEL_DIR}/model_dirac_JBand.cpp
${MODEL_DIR}/model_dirac_JBandB.cpp
${MODEL_DIR}/model_dirac_JDiag.cpp
${MODEL_DIR}/model_dirac_JSparse.cpp
${MODEL_DIR}/model_dirac_JSparseB.cpp
${MODEL_DIR}/model_dirac_Jrz.cpp
${MODEL_DIR}/model_dirac_Jv.cpp
${MODEL_DIR}/model_dirac_JvB.cpp
${MODEL_DIR}/model_dirac_Jy.cpp
${MODEL_DIR}/model_dirac_Jz.cpp
${MODEL_DIR}/model_dirac_dJrzdsigma.cpp
${MODEL_DIR}/model_dirac_dJrzdz.cpp
${MODEL_DIR}/model_dirac_dJydsigma.cpp
${MODEL_DIR}/model_dirac_dJydy.cpp
${MODEL_DIR}/model_dirac_dJzdsigma.cpp
${MODEL_DIR}/model_dirac_dJzdz.cpp
${MODEL_DIR}/model_dirac_deltaqB.cpp
${MODEL_DIR}/model_dirac_deltasx.cpp
${MODEL_DIR}/model_dirac_deltax.cpp
${MODEL_DIR}/model_dirac_deltaxB.cpp
${MODEL_DIR}/model_dirac_drzdp.cpp
${MODEL_DIR}/model_dirac_drzdx.cpp
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
${MODEL_DIR}/model_dirac_rz.cpp
${MODEL_DIR}/model_dirac_sigma_y.cpp
${MODEL_DIR}/model_dirac_sigma_z.cpp
${MODEL_DIR}/model_dirac_srz.cpp
${MODEL_DIR}/model_dirac_stau.cpp
${MODEL_DIR}/model_dirac_sx0.cpp
${MODEL_DIR}/model_dirac_sxdot.cpp
${MODEL_DIR}/model_dirac_sz.cpp
${MODEL_DIR}/model_dirac_w.cpp
${MODEL_DIR}/model_dirac_x0.cpp
${MODEL_DIR}/model_dirac_xBdot.cpp
${MODEL_DIR}/model_dirac_xdot.cpp
${MODEL_DIR}/model_dirac_y.cpp
${MODEL_DIR}/model_dirac_z.cpp

${MODEL_DIR}/wrapfunctions.cpp
)
    
add_library(${PROJECT_NAME} ${SRC_LIST_LIB})

if(APPLE)
    set(SHARED_OBJECT_EXTENSION .dylib)
else()
    set(SHARED_OBJECT_EXTENSION .so)
endif()
    
target_link_libraries(
    ${PROJECT_NAME}
    Upstream::amici
)

set(SRC_LIST_EXE main.cpp)

add_executable(simulate_${PROJECT_NAME} ${SRC_LIST_EXE})
    
target_link_libraries(simulate_${PROJECT_NAME} Upstream::amici ${PROJECT_NAME})

if($ENV{ENABLE_GCOV_COVERAGE})
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 --coverage")
endif()
