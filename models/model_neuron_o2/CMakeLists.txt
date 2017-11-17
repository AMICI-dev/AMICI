project(model_neuron_o2)

cmake_minimum_required(VERSION 2.8)

include(CheckCXXCompilerFlag)
set(MY_CXX_FLAGS -std=c++11 -Wall -Wno-unused-function -Wno-unused-variable -Wno-unused-but-set-variable)
foreach(FLAG ${MY_CXX_FLAGS})
    unset(CUR_FLAG_SUPPORTED CACHE)
    CHECK_CXX_COMPILER_FLAG(${FLAG} CUR_FLAG_SUPPORTED)
    if(${CUR_FLAG_SUPPORTED})
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAG}")
    endif()
endforeach(FLAG)

find_package(Amici HINTS ${CMAKE_CURRENT_LIST_DIR}/../../build)
get_property(AmiciConfigIncludes TARGET Upstream::amici PROPERTY INCLUDE_DIRECTORIES)

set(MODEL_DIR "${Amici_DIR}/../models/model_neuron_o2")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories(${AmiciConfigIncludes})

set(SRC_LIST_LIB ${MODEL_DIR}/model_neuron_o2_J.cpp
${MODEL_DIR}/model_neuron_o2_JB.cpp
${MODEL_DIR}/model_neuron_o2_JDiag.cpp
${MODEL_DIR}/model_neuron_o2_JSparse.cpp
${MODEL_DIR}/model_neuron_o2_JSparseB.cpp
${MODEL_DIR}/model_neuron_o2_Jrz.cpp
${MODEL_DIR}/model_neuron_o2_Jv.cpp
${MODEL_DIR}/model_neuron_o2_JvB.cpp
${MODEL_DIR}/model_neuron_o2_Jy.cpp
${MODEL_DIR}/model_neuron_o2_Jz.cpp
${MODEL_DIR}/model_neuron_o2_dJrzdsigma.cpp
${MODEL_DIR}/model_neuron_o2_dJrzdz.cpp
${MODEL_DIR}/model_neuron_o2_dJydsigma.cpp
${MODEL_DIR}/model_neuron_o2_dJydy.cpp
${MODEL_DIR}/model_neuron_o2_dJzdsigma.cpp
${MODEL_DIR}/model_neuron_o2_dJzdz.cpp
${MODEL_DIR}/model_neuron_o2_deltaqB.cpp
${MODEL_DIR}/model_neuron_o2_deltasx.cpp
${MODEL_DIR}/model_neuron_o2_deltax.cpp
${MODEL_DIR}/model_neuron_o2_deltaxB.cpp
${MODEL_DIR}/model_neuron_o2_drzdp.cpp
${MODEL_DIR}/model_neuron_o2_drzdx.cpp
${MODEL_DIR}/model_neuron_o2_dsigma_ydp.cpp
${MODEL_DIR}/model_neuron_o2_dsigma_zdp.cpp
${MODEL_DIR}/model_neuron_o2_dwdp.cpp
${MODEL_DIR}/model_neuron_o2_dwdx.cpp
${MODEL_DIR}/model_neuron_o2_dxdotdp.cpp
${MODEL_DIR}/model_neuron_o2_dydp.cpp
${MODEL_DIR}/model_neuron_o2_dydx.cpp
${MODEL_DIR}/model_neuron_o2_dzdp.cpp
${MODEL_DIR}/model_neuron_o2_dzdx.cpp
${MODEL_DIR}/model_neuron_o2_qBdot.cpp
${MODEL_DIR}/model_neuron_o2_root.cpp
${MODEL_DIR}/model_neuron_o2_rz.cpp
${MODEL_DIR}/model_neuron_o2_sigma_y.cpp
${MODEL_DIR}/model_neuron_o2_sigma_z.cpp
${MODEL_DIR}/model_neuron_o2_srz.cpp
${MODEL_DIR}/model_neuron_o2_stau.cpp
${MODEL_DIR}/model_neuron_o2_sx0.cpp
${MODEL_DIR}/model_neuron_o2_sxdot.cpp
${MODEL_DIR}/model_neuron_o2_sz.cpp
${MODEL_DIR}/model_neuron_o2_w.cpp
${MODEL_DIR}/model_neuron_o2_x0.cpp
${MODEL_DIR}/model_neuron_o2_xBdot.cpp
${MODEL_DIR}/model_neuron_o2_xdot.cpp
${MODEL_DIR}/model_neuron_o2_y.cpp
${MODEL_DIR}/model_neuron_o2_z.cpp

${MODEL_DIR}/wrapfunctions.cpp
)
    
add_library(${PROJECT_NAME} ${SRC_LIST_LIB})
    
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
