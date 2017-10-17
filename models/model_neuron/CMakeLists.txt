project(model_neuron)

cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wno-unused-function -Wno-unused-variable -Wno-unused-but-set-variable")

find_package(Amici)
get_property(AmiciConfigIncludes TARGET Upstream::amici PROPERTY INCLUDE_DIRECTORIES)

set(MODEL_DIR "${Amici_DIR}/../models/model_neuron")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories(${AmiciConfigIncludes})

set(SRC_LIST_LIB ${MODEL_DIR}/model_neuron_J.cpp
${MODEL_DIR}/model_neuron_JB.cpp
${MODEL_DIR}/model_neuron_JBand.cpp
${MODEL_DIR}/model_neuron_JBandB.cpp
${MODEL_DIR}/model_neuron_JDiag.cpp
${MODEL_DIR}/model_neuron_JSparse.cpp
${MODEL_DIR}/model_neuron_JSparseB.cpp
${MODEL_DIR}/model_neuron_Jrz.cpp
${MODEL_DIR}/model_neuron_Jv.cpp
${MODEL_DIR}/model_neuron_JvB.cpp
${MODEL_DIR}/model_neuron_Jy.cpp
${MODEL_DIR}/model_neuron_Jz.cpp
${MODEL_DIR}/model_neuron_dJrzdsigma.cpp
${MODEL_DIR}/model_neuron_dJrzdz.cpp
${MODEL_DIR}/model_neuron_dJydsigma.cpp
${MODEL_DIR}/model_neuron_dJydy.cpp
${MODEL_DIR}/model_neuron_dJzdsigma.cpp
${MODEL_DIR}/model_neuron_dJzdz.cpp
${MODEL_DIR}/model_neuron_deltaqB.cpp
${MODEL_DIR}/model_neuron_deltasx.cpp
${MODEL_DIR}/model_neuron_deltax.cpp
${MODEL_DIR}/model_neuron_deltaxB.cpp
${MODEL_DIR}/model_neuron_drzdp.cpp
${MODEL_DIR}/model_neuron_drzdx.cpp
${MODEL_DIR}/model_neuron_dsigma_ydp.cpp
${MODEL_DIR}/model_neuron_dsigma_zdp.cpp
${MODEL_DIR}/model_neuron_dwdp.cpp
${MODEL_DIR}/model_neuron_dwdx.cpp
${MODEL_DIR}/model_neuron_dxdotdp.cpp
${MODEL_DIR}/model_neuron_dydp.cpp
${MODEL_DIR}/model_neuron_dydx.cpp
${MODEL_DIR}/model_neuron_dzdp.cpp
${MODEL_DIR}/model_neuron_dzdx.cpp
${MODEL_DIR}/model_neuron_qBdot.cpp
${MODEL_DIR}/model_neuron_root.cpp
${MODEL_DIR}/model_neuron_rz.cpp
${MODEL_DIR}/model_neuron_sigma_y.cpp
${MODEL_DIR}/model_neuron_sigma_z.cpp
${MODEL_DIR}/model_neuron_srz.cpp
${MODEL_DIR}/model_neuron_stau.cpp
${MODEL_DIR}/model_neuron_sx0.cpp
${MODEL_DIR}/model_neuron_sxdot.cpp
${MODEL_DIR}/model_neuron_sz.cpp
${MODEL_DIR}/model_neuron_w.cpp
${MODEL_DIR}/model_neuron_x0.cpp
${MODEL_DIR}/model_neuron_xBdot.cpp
${MODEL_DIR}/model_neuron_xdot.cpp
${MODEL_DIR}/model_neuron_y.cpp
${MODEL_DIR}/model_neuron_z.cpp

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
