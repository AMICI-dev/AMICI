project(model_dirac)

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

set(MODEL_DIR "${Amici_DIR}/../models/model_dirac")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories(${AmiciConfigIncludes})

set(SRC_LIST_LIB ${MODEL_DIR}/model_dirac_J.cpp
${MODEL_DIR}/model_dirac_JB.cpp
${MODEL_DIR}/model_dirac_JDiag.cpp
${MODEL_DIR}/model_dirac_JSparse.cpp
${MODEL_DIR}/model_dirac_JSparseB.cpp
${MODEL_DIR}/model_dirac_Jv.cpp
${MODEL_DIR}/model_dirac_JvB.cpp
${MODEL_DIR}/model_dirac_Jy.cpp
${MODEL_DIR}/model_dirac_dJydsigma.cpp
${MODEL_DIR}/model_dirac_dJydy.cpp
${MODEL_DIR}/model_dirac_deltasx.cpp
${MODEL_DIR}/model_dirac_deltax.cpp
${MODEL_DIR}/model_dirac_dxdotdp.cpp
${MODEL_DIR}/model_dirac_dydx.cpp
${MODEL_DIR}/model_dirac_qBdot.cpp
${MODEL_DIR}/model_dirac_root.cpp
${MODEL_DIR}/model_dirac_sigma_y.cpp
${MODEL_DIR}/model_dirac_stau.cpp
${MODEL_DIR}/model_dirac_sxdot.cpp
${MODEL_DIR}/model_dirac_xBdot.cpp
${MODEL_DIR}/model_dirac_xdot.cpp
${MODEL_DIR}/model_dirac_y.cpp

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
