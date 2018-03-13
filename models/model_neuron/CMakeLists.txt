project(model_neuron)

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

set(MODEL_DIR "${Amici_DIR}/../models/model_neuron")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories(${AmiciConfigIncludes})

set(SRC_LIST_LIB ${MODEL_DIR}/model_neuron_J.cpp
${MODEL_DIR}/model_neuron_JB.cpp
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
${MODEL_DIR}/model_neuron_drzdx.cpp
${MODEL_DIR}/model_neuron_dxdotdp.cpp
${MODEL_DIR}/model_neuron_dydx.cpp
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
${MODEL_DIR}/model_neuron_x0.cpp
${MODEL_DIR}/model_neuron_xBdot.cpp
${MODEL_DIR}/model_neuron_xdot.cpp
${MODEL_DIR}/model_neuron_y.cpp
${MODEL_DIR}/model_neuron_z.cpp

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

## SWIG

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

FIND_PACKAGE(PythonLibs REQUIRED)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

SET(CMAKE_SWIG_FLAGS "")
SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_LIBRARY(${PROJECT_NAME} MODULE LANGUAGE python SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}.i)
SWIG_LINK_LIBRARIES(${PROJECT_NAME}
	Upstream::amici
	${PYTHON_LIBRARIES} 
	${CMAKE_CURRENT_BINARY_DIR}/lib${PROJECT_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX})


SET(SETUP_PY_IN ${AMICI_DIR}/model_setup.template.py)
SET(SETUP_PY_OUT ${CMAKE_CURRENT_BINARY_DIR}/setup.py)

add_custom_target(install-python
        DEPENDS _${PROJECT_NAME}
        COMMAND python ${SETUP_PY_OUT} install)



