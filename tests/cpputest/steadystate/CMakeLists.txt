project(model_steadystate_test)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wno-unused-function")
add_definitions(-DAMICI_WITHOUT_MATLAB)

set(AMICI_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../")
set(MODEL_DIR "${AMICI_DIR}/models/model_steadystate")
set(SUITESPARSE_DIR "${AMICI_DIR}/SuiteSparse/")
set(SUITESPARSE_LIB_DIR "${AMICI_DIR}/SuiteSparse/lib")
set(SUNDIALS_LIB_DIR "${AMICI_DIR}/sundials/build/lib")

set(CPPUTEST_DIR "${AMICI_DIR}/ThirdParty/cpputest-3.8/")
set(CPPUTEST_INCLUDE_DIR "${CPPUTEST_DIR}/include")
set(CPPUTEST_LIBRARY "${CPPUTEST_DIR}/cpputest_build/lib/libCppUTest.a" "${CPPUTEST_DIR}/cpputest_build/lib/libCppUTestExt.a")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -include ${CPPUTEST_INCLUDE_DIR}/CppUTest/MemoryLeakDetectorMallocMacros.h")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -include ${CPPUTEST_INCLUDE_DIR}/CppUTest/MemoryLeakDetectorNewMacros.h")

find_package(HDF5 COMPONENTS C HL REQUIRED)

set(SRC_LIST
    main.cpp
    tests1.cpp
    ../testfunctions.cpp
)

include_directories("${AMICI_DIR}")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/..")
include_directories("${MODEL_DIR}")
include_directories("${HDF5_INCLUDE_DIRS}")
include_directories("${AMICI_DIR}/sundials/include")
include_directories("${SUITESPARSE_DIR}/include")
include_directories(${CPPUTEST_INCLUDE_DIR})

add_executable(${PROJECT_NAME} ${SRC_LIST})

if(APPLE)
    set(SHARED_OBJECT_EXTENSION .dylib)
else()
    set(SHARED_OBJECT_EXTENSION .so)
endif()

target_link_libraries(${PROJECT_NAME}
${AMICI_DIR}/models/model_steadystate/build/libmodel_steadystate.a
${CPPUTEST_LIBRARY}
"${SUNDIALS_LIB_DIR}/libsundials_nvecserial${SHARED_OBJECT_EXTENSION}"
"${SUNDIALS_LIB_DIR}/libsundials_cvodes${SHARED_OBJECT_EXTENSION}"
"${SUITESPARSE_LIB_DIR}/libcolamd${SHARED_OBJECT_EXTENSION}"
"${SUITESPARSE_LIB_DIR}/libklu${SHARED_OBJECT_EXTENSION}"
"${SUITESPARSE_LIB_DIR}/libbtf${SHARED_OBJECT_EXTENSION}"
"${SUITESPARSE_LIB_DIR}/libamd${SHARED_OBJECT_EXTENSION}"
"${SUITESPARSE_LIB_DIR}/libsuitesparseconfig${SHARED_OBJECT_EXTENSION}"

"${HDF5_HL_LIBRARIES}"
"${HDF5_C_LIBRARIES}"
"-lpthread -ldl -lz"
"-lm"
)
