project(amiciIntegrationTests)

cmake_minimum_required(VERSION 2.8)

set(CMAKE_BUILD_TYPE Debug)

find_package(Amici HINTS ${CMAKE_CURRENT_LIST_DIR}/../../build)
get_property(AmiciConfigIncludes TARGET Upstream::amici PROPERTY INCLUDE_DIRECTORIES)
set(AMICI_DIR ${Amici_DIR}/..)

set(CPPUTEST_DIR "${AMICI_DIR}/ThirdParty/cpputest-3.8/" CACHE PATH "cppUTest root directory")
set(CPPUTEST_INCLUDE_DIR "${CPPUTEST_DIR}/include" CACHE PATH "cppUTest include directory")
set(CPPUTEST_LIBRARY "${CPPUTEST_DIR}/lib/libCppUTest.a" "${CPPUTEST_DIR}/lib/libCppUTestExt.a"  CACHE PATH "cppUTest libraries")

set(CMAKE_CXX_FLAGS_OLD "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -include sstream -include ${CPPUTEST_INCLUDE_DIR}/CppUTest/MemoryLeakDetectorNewMacros.h")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories(${AmiciConfigIncludes})
include_directories(${CPPUTEST_INCLUDE_DIR})

if(APPLE)
    set(SHARED_OBJECT_EXTENSION .dylib)
else()
    set(SHARED_OBJECT_EXTENSION .so)
endif()

enable_testing()

add_subdirectory(unittests)
add_subdirectory(dirac)
add_subdirectory(events)
add_subdirectory(nested_events)
add_subdirectory(steadystate)
add_subdirectory(jakstat_adjoint)
add_subdirectory(neuron)
add_subdirectory(jakstat_adjoint_o2)
add_subdirectory(neuron_o2)

if($ENV{ENABLE_GCOV_COVERAGE})
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 --coverage")
endif()
