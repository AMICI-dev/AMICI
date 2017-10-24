project(model_steadystate)

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

find_package(Amici)
get_property(AmiciConfigIncludes TARGET Upstream::amici PROPERTY INCLUDE_DIRECTORIES)

set(MODEL_DIR "${Amici_DIR}/../models/model_steadystate")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories(${AmiciConfigIncludes})

set(SRC_LIST_LIB ${MODEL_DIR}/model_steadystate_xdot.cpp
${MODEL_DIR}/model_steadystate_w.cpp
${MODEL_DIR}/model_steadystate_dwdx.cpp
${MODEL_DIR}/model_steadystate_J.cpp
${MODEL_DIR}/model_steadystate_x0.cpp
${MODEL_DIR}/model_steadystate_Jv.cpp
${MODEL_DIR}/model_steadystate_JBand.cpp
${MODEL_DIR}/model_steadystate_JSparse.cpp
${MODEL_DIR}/model_steadystate_JDiag.cpp
${MODEL_DIR}/model_steadystate_y.cpp
${MODEL_DIR}/model_steadystate_z.cpp
${MODEL_DIR}/model_steadystate_rz.cpp
${MODEL_DIR}/model_steadystate_deltax.cpp
${MODEL_DIR}/model_steadystate_root.cpp
${MODEL_DIR}/model_steadystate_Jy.cpp
${MODEL_DIR}/model_steadystate_Jz.cpp
${MODEL_DIR}/model_steadystate_Jrz.cpp
${MODEL_DIR}/model_steadystate_sigma_y.cpp
${MODEL_DIR}/model_steadystate_sigma_z.cpp
${MODEL_DIR}/model_steadystate_sxdot.cpp
${MODEL_DIR}/model_steadystate_sx0.cpp
${MODEL_DIR}/model_steadystate_sz.cpp
${MODEL_DIR}/model_steadystate_deltasx.cpp
${MODEL_DIR}/model_steadystate_stau.cpp
${MODEL_DIR}/model_steadystate_srz.cpp
${MODEL_DIR}/model_steadystate_dJydy.cpp
${MODEL_DIR}/model_steadystate_dJydsigma.cpp
${MODEL_DIR}/model_steadystate_dJzdz.cpp
${MODEL_DIR}/model_steadystate_dJzdsigma.cpp
${MODEL_DIR}/model_steadystate_dJrzdz.cpp
${MODEL_DIR}/model_steadystate_dJrzdsigma.cpp
${MODEL_DIR}/model_steadystate_dwdp.cpp
${MODEL_DIR}/model_steadystate_dxdotdp.cpp
${MODEL_DIR}/model_steadystate_dydp.cpp
${MODEL_DIR}/model_steadystate_dsigma_ydp.cpp
${MODEL_DIR}/model_steadystate_dsigma_zdp.cpp
${MODEL_DIR}/model_steadystate_dydx.cpp
${MODEL_DIR}/model_steadystate_dzdx.cpp
${MODEL_DIR}/model_steadystate_dzdp.cpp
${MODEL_DIR}/model_steadystate_drzdx.cpp
${MODEL_DIR}/model_steadystate_drzdp.cpp
${MODEL_DIR}/model_steadystate_xBdot.cpp
${MODEL_DIR}/model_steadystate_qBdot.cpp
${MODEL_DIR}/model_steadystate_JB.cpp
${MODEL_DIR}/model_steadystate_JvB.cpp
${MODEL_DIR}/model_steadystate_JBandB.cpp
${MODEL_DIR}/model_steadystate_JSparseB.cpp
${MODEL_DIR}/model_steadystate_deltaxB.cpp
${MODEL_DIR}/model_steadystate_deltaqB.cpp

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
