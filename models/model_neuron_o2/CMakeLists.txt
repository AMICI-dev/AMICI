cmake_minimum_required(VERSION 3.15)

project(model_neuron_o2)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(CheckCXXCompilerFlag)
set(MY_CXX_FLAGS -Wall -Wno-unused-function -Wno-unused-variable)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  list(APPEND MY_CXX_FLAGS -Wno-unused-but-set-variable)
endif()
foreach(FLAG ${MY_CXX_FLAGS})
  unset(CUR_FLAG_SUPPORTED CACHE)
  check_cxx_compiler_flag(-Werror ${FLAG} CUR_FLAG_SUPPORTED)
  if(${CUR_FLAG_SUPPORTED})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAG}")
  endif()
endforeach(FLAG)

find_package(Amici  REQUIRED HINTS
             ${CMAKE_CURRENT_LIST_DIR}/../../build)
message(STATUS "Found AMICI ${Amici_DIR}")

set(MODEL_DIR ${CMAKE_CURRENT_LIST_DIR})

set(SRC_LIST_LIB ${MODEL_DIR}/model_neuron_o2_JSparse.cpp
${MODEL_DIR}/model_neuron_o2_Jrz.cpp
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
${MODEL_DIR}/model_neuron_o2_drzdx.cpp
${MODEL_DIR}/model_neuron_o2_dwdx.cpp
${MODEL_DIR}/model_neuron_o2_dxdotdp.cpp
${MODEL_DIR}/model_neuron_o2_dydx.cpp
${MODEL_DIR}/model_neuron_o2_dzdx.cpp
${MODEL_DIR}/model_neuron_o2_root.cpp
${MODEL_DIR}/model_neuron_o2_rz.cpp
${MODEL_DIR}/model_neuron_o2_sigmay.cpp
${MODEL_DIR}/model_neuron_o2_sigmaz.cpp
${MODEL_DIR}/model_neuron_o2_srz.cpp
${MODEL_DIR}/model_neuron_o2_stau.cpp
${MODEL_DIR}/model_neuron_o2_sx0.cpp
${MODEL_DIR}/model_neuron_o2_sz.cpp
${MODEL_DIR}/model_neuron_o2_w.cpp
${MODEL_DIR}/model_neuron_o2_x0.cpp
${MODEL_DIR}/model_neuron_o2_xdot.cpp
${MODEL_DIR}/model_neuron_o2_y.cpp
${MODEL_DIR}/model_neuron_o2_z.cpp
 ${MODEL_DIR}/wrapfunctions.cpp)

add_library(${PROJECT_NAME} ${SRC_LIST_LIB})
add_library(model ALIAS ${PROJECT_NAME})

# This option can be helpful when using the Intel compiler and compilation of
# wrapfunctions.cpp fails due to insufficient memory.
option(ENABLE_WRAPFUNCTIONS_OPTIMIZATIONS
       "Enable compiler optimizations for wrapfunctions.cpp?" ON)
if(NOT ENABLE_WRAPFUNCTIONS_OPTIMIZATIONS)
  set_source_files_properties(wrapfunctions.cpp PROPERTIES COMPILE_FLAGS -O0)
endif()

# Some special functions require boost
# TODO: set some flag during code generation whether the given model requires
# boost. for now, try to find it, add include directories and link against it.
# let the compiler/linker error if it is required but not found
find_package(Boost COMPONENTS math)

target_include_directories(${PROJECT_NAME}
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
    PRIVATE $<$<BOOL:${Boost_FOUND}>:Boost::math>
)

target_link_libraries(${PROJECT_NAME} PUBLIC Upstream::amici)

set(SRC_LIST_EXE main.cpp)

add_executable(simulate_${PROJECT_NAME} ${SRC_LIST_EXE})

target_link_libraries(simulate_${PROJECT_NAME} ${PROJECT_NAME})

if($ENV{ENABLE_GCOV_COVERAGE})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 --coverage")
endif()

# SWIG
option(ENABLE_SWIG "Build swig/python library?" ON)
if(ENABLE_SWIG)
  if(NOT (${CMAKE_VERSION} VERSION_LESS 3.8))
    add_subdirectory(swig)
  else()
    message(WARNING "Unable to build SWIG interface, upgrade CMake to >=3.8.")
  endif()
endif()

# <Export cmake configuration>
include(GNUInstallDirs)
install(
  TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}Targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
export(
  EXPORT ${PROJECT_NAME}Targets
  FILE ${PROJECT_NAME}Config.cmake
  NAMESPACE Upstream::)
# </Export cmake configuration>
