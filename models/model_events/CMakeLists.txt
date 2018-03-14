project(model_events)

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

set(MODEL_DIR "${Amici_DIR}/../models/model_events")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories(${AmiciConfigIncludes})

set(SRC_LIST_LIB ${MODEL_DIR}/model_events_J.cpp
${MODEL_DIR}/model_events_JB.cpp
${MODEL_DIR}/model_events_JDiag.cpp
${MODEL_DIR}/model_events_JSparse.cpp
${MODEL_DIR}/model_events_JSparseB.cpp
${MODEL_DIR}/model_events_Jrz.cpp
${MODEL_DIR}/model_events_Jv.cpp
${MODEL_DIR}/model_events_JvB.cpp
${MODEL_DIR}/model_events_Jy.cpp
${MODEL_DIR}/model_events_Jz.cpp
${MODEL_DIR}/model_events_dJrzdsigma.cpp
${MODEL_DIR}/model_events_dJrzdz.cpp
${MODEL_DIR}/model_events_dJydsigma.cpp
${MODEL_DIR}/model_events_dJydy.cpp
${MODEL_DIR}/model_events_dJzdsigma.cpp
${MODEL_DIR}/model_events_dJzdz.cpp
${MODEL_DIR}/model_events_deltasx.cpp
${MODEL_DIR}/model_events_drzdx.cpp
${MODEL_DIR}/model_events_dxdotdp.cpp
${MODEL_DIR}/model_events_dydp.cpp
${MODEL_DIR}/model_events_dydx.cpp
${MODEL_DIR}/model_events_dzdx.cpp
${MODEL_DIR}/model_events_qBdot.cpp
${MODEL_DIR}/model_events_root.cpp
${MODEL_DIR}/model_events_rz.cpp
${MODEL_DIR}/model_events_sigma_y.cpp
${MODEL_DIR}/model_events_sigma_z.cpp
${MODEL_DIR}/model_events_srz.cpp
${MODEL_DIR}/model_events_stau.cpp
${MODEL_DIR}/model_events_sxdot.cpp
${MODEL_DIR}/model_events_sz.cpp
${MODEL_DIR}/model_events_x0.cpp
${MODEL_DIR}/model_events_xBdot.cpp
${MODEL_DIR}/model_events_xdot.cpp
${MODEL_DIR}/model_events_y.cpp
${MODEL_DIR}/model_events_z.cpp

${MODEL_DIR}/wrapfunctions.cpp
)
    
add_library(${PROJECT_NAME} ${SRC_LIST_LIB})
    
target_link_libraries(
    ${PROJECT_NAME}
    Upstream::amici
)

set(SRC_LIST_EXE main.cpp)

add_executable(simulate_${PROJECT_NAME} ${SRC_LIST_EXE})
    
target_link_libraries(simulate_${PROJECT_NAME} ${PROJECT_NAME} Upstream::amici)

if($ENV{ENABLE_GCOV_COVERAGE})
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 --coverage")
endif()
