project(model_steadystate)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(CheckCXXCompilerFlag)
set(MY_CXX_FLAGS -Wall -Wno-unused-function -Wno-unused-variable -Wno-unused-but-set-variable)
foreach(FLAG ${MY_CXX_FLAGS})
    unset(CUR_FLAG_SUPPORTED CACHE)
    CHECK_CXX_COMPILER_FLAG(${FLAG} CUR_FLAG_SUPPORTED)
    if(${CUR_FLAG_SUPPORTED})
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAG}")
    endif()
endforeach(FLAG)

find_package(Amici HINTS ${CMAKE_CURRENT_LIST_DIR}/../../build)

set(MODEL_DIR ${CMAKE_CURRENT_LIST_DIR})

set(SRC_LIST_LIB ${MODEL_DIR}/model_steadystate_J.cpp
${MODEL_DIR}/model_steadystate_JB.cpp
${MODEL_DIR}/model_steadystate_JDiag.cpp
${MODEL_DIR}/model_steadystate_JSparse.cpp
${MODEL_DIR}/model_steadystate_JSparseB.cpp
${MODEL_DIR}/model_steadystate_Jv.cpp
${MODEL_DIR}/model_steadystate_JvB.cpp
${MODEL_DIR}/model_steadystate_Jy.cpp
${MODEL_DIR}/model_steadystate_dJydsigma.cpp
${MODEL_DIR}/model_steadystate_dJydy.cpp
${MODEL_DIR}/model_steadystate_dwdp.cpp
${MODEL_DIR}/model_steadystate_dwdx.cpp
${MODEL_DIR}/model_steadystate_dxdotdp.cpp
${MODEL_DIR}/model_steadystate_dydx.cpp
${MODEL_DIR}/model_steadystate_qBdot.cpp
${MODEL_DIR}/model_steadystate_sigma_y.cpp
${MODEL_DIR}/model_steadystate_sxdot.cpp
${MODEL_DIR}/model_steadystate_w.cpp
${MODEL_DIR}/model_steadystate_x0.cpp
${MODEL_DIR}/model_steadystate_xBdot.cpp
${MODEL_DIR}/model_steadystate_xdot.cpp
${MODEL_DIR}/model_steadystate_y.cpp

${MODEL_DIR}/wrapfunctions.cpp
)

add_library(${PROJECT_NAME} ${SRC_LIST_LIB})
add_library(model ALIAS ${PROJECT_NAME})

target_include_directories(${PROJECT_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
    
target_link_libraries(${PROJECT_NAME}
    PUBLIC Upstream::amici
)

set(SRC_LIST_EXE main.cpp)

add_executable(simulate_${PROJECT_NAME} ${SRC_LIST_EXE})
    
target_link_libraries(simulate_${PROJECT_NAME} ${PROJECT_NAME})

if($ENV{ENABLE_GCOV_COVERAGE})
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 --coverage")
endif()

## SWIG

add_subdirectory(swig)

# <Export cmake configuration>
include(GNUInstallDirs)
install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}Targets
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
export(EXPORT ${PROJECT_NAME}Targets FILE ${PROJECT_NAME}Config.cmake
    NAMESPACE Upstream::
    )
# </Export cmake configuration>

