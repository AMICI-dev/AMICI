project(model_nested_events)

cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wno-unused-function -Wno-unused-variable -Wno-unused-but-set-variable")

find_package(Amici HINTS ${CMAKE_CURRENT_LIST_DIR}/../../build)
get_property(AmiciConfigIncludes TARGET Upstream::amici PROPERTY INCLUDE_DIRECTORIES)

set(MODEL_DIR "${Amici_DIR}/../models/model_nested_events")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories(${AmiciConfigIncludes})

set(SRC_LIST_LIB ${MODEL_DIR}/model_nested_events_J.cpp
${MODEL_DIR}/model_nested_events_JB.cpp
${MODEL_DIR}/model_nested_events_JBand.cpp
${MODEL_DIR}/model_nested_events_JBandB.cpp
${MODEL_DIR}/model_nested_events_JDiag.cpp
${MODEL_DIR}/model_nested_events_JSparse.cpp
${MODEL_DIR}/model_nested_events_JSparseB.cpp
${MODEL_DIR}/model_nested_events_Jrz.cpp
${MODEL_DIR}/model_nested_events_Jv.cpp
${MODEL_DIR}/model_nested_events_JvB.cpp
${MODEL_DIR}/model_nested_events_Jy.cpp
${MODEL_DIR}/model_nested_events_Jz.cpp
${MODEL_DIR}/model_nested_events_dJrzdsigma.cpp
${MODEL_DIR}/model_nested_events_dJrzdz.cpp
${MODEL_DIR}/model_nested_events_dJydsigma.cpp
${MODEL_DIR}/model_nested_events_dJydy.cpp
${MODEL_DIR}/model_nested_events_dJzdsigma.cpp
${MODEL_DIR}/model_nested_events_dJzdz.cpp
${MODEL_DIR}/model_nested_events_deltaqB.cpp
${MODEL_DIR}/model_nested_events_deltasx.cpp
${MODEL_DIR}/model_nested_events_deltax.cpp
${MODEL_DIR}/model_nested_events_deltaxB.cpp
${MODEL_DIR}/model_nested_events_drzdp.cpp
${MODEL_DIR}/model_nested_events_drzdx.cpp
${MODEL_DIR}/model_nested_events_dsigma_ydp.cpp
${MODEL_DIR}/model_nested_events_dsigma_zdp.cpp
${MODEL_DIR}/model_nested_events_dwdp.cpp
${MODEL_DIR}/model_nested_events_dwdx.cpp
${MODEL_DIR}/model_nested_events_dxdotdp.cpp
${MODEL_DIR}/model_nested_events_dydp.cpp
${MODEL_DIR}/model_nested_events_dydx.cpp
${MODEL_DIR}/model_nested_events_dzdp.cpp
${MODEL_DIR}/model_nested_events_dzdx.cpp
${MODEL_DIR}/model_nested_events_qBdot.cpp
${MODEL_DIR}/model_nested_events_root.cpp
${MODEL_DIR}/model_nested_events_rz.cpp
${MODEL_DIR}/model_nested_events_sigma_y.cpp
${MODEL_DIR}/model_nested_events_sigma_z.cpp
${MODEL_DIR}/model_nested_events_srz.cpp
${MODEL_DIR}/model_nested_events_stau.cpp
${MODEL_DIR}/model_nested_events_sx0.cpp
${MODEL_DIR}/model_nested_events_sxdot.cpp
${MODEL_DIR}/model_nested_events_sz.cpp
${MODEL_DIR}/model_nested_events_w.cpp
${MODEL_DIR}/model_nested_events_x0.cpp
${MODEL_DIR}/model_nested_events_xBdot.cpp
${MODEL_DIR}/model_nested_events_xdot.cpp
${MODEL_DIR}/model_nested_events_y.cpp
${MODEL_DIR}/model_nested_events_z.cpp

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
