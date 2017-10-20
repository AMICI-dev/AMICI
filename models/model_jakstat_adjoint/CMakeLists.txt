project(model_jakstat_adjoint)

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

set(MODEL_DIR "${Amici_DIR}/../models/model_jakstat_adjoint")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories(${AmiciConfigIncludes})

set(SRC_LIST_LIB ${MODEL_DIR}/model_jakstat_adjoint_J.cpp
${MODEL_DIR}/model_jakstat_adjoint_JB.cpp
${MODEL_DIR}/model_jakstat_adjoint_JBand.cpp
${MODEL_DIR}/model_jakstat_adjoint_JBandB.cpp
${MODEL_DIR}/model_jakstat_adjoint_JDiag.cpp
${MODEL_DIR}/model_jakstat_adjoint_JSparse.cpp
${MODEL_DIR}/model_jakstat_adjoint_JSparseB.cpp
${MODEL_DIR}/model_jakstat_adjoint_Jrz.cpp
${MODEL_DIR}/model_jakstat_adjoint_Jv.cpp
${MODEL_DIR}/model_jakstat_adjoint_JvB.cpp
${MODEL_DIR}/model_jakstat_adjoint_Jy.cpp
${MODEL_DIR}/model_jakstat_adjoint_Jz.cpp
${MODEL_DIR}/model_jakstat_adjoint_dJrzdsigma.cpp
${MODEL_DIR}/model_jakstat_adjoint_dJrzdz.cpp
${MODEL_DIR}/model_jakstat_adjoint_dJydsigma.cpp
${MODEL_DIR}/model_jakstat_adjoint_dJydy.cpp
${MODEL_DIR}/model_jakstat_adjoint_dJzdsigma.cpp
${MODEL_DIR}/model_jakstat_adjoint_dJzdz.cpp
${MODEL_DIR}/model_jakstat_adjoint_deltaqB.cpp
${MODEL_DIR}/model_jakstat_adjoint_deltasx.cpp
${MODEL_DIR}/model_jakstat_adjoint_deltax.cpp
${MODEL_DIR}/model_jakstat_adjoint_deltaxB.cpp
${MODEL_DIR}/model_jakstat_adjoint_drzdp.cpp
${MODEL_DIR}/model_jakstat_adjoint_drzdx.cpp
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
${MODEL_DIR}/model_jakstat_adjoint_rz.cpp
${MODEL_DIR}/model_jakstat_adjoint_sigma_y.cpp
${MODEL_DIR}/model_jakstat_adjoint_sigma_z.cpp
${MODEL_DIR}/model_jakstat_adjoint_srz.cpp
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
