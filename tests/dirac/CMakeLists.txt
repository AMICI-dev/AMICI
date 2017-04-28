project(model_dirac_test)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wno-unused-function")
add_definitions(-DAMICI_WITHOUT_MATLAB)

set(AMICI_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../")
set(MODEL_DIR "${AMICI_DIR}/models/model_dirac")
set(SUITESPARSE_DIR "${AMICI_DIR}/SuiteSparse/")
set(SUITESPARSE_LIB_DIR "${AMICI_DIR}/SuiteSparse/")
set(SUNDIALS_LIB_DIR "${AMICI_DIR}/build/sundials/lib")

set(CPPUTEST_DIR "/home/dweindl/src/cpputest/")
set(CPPUTEST_INCLUDE_DIR "${CPPUTEST_DIR}/include")
set(CPPUTEST_LIBRARY "${CPPUTEST_DIR}/build/src/CppUTest/libCppUTest.a" "${CPPUTEST_DIR}/build/src/CppUTestExt/libCppUTestExt.a")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -include ${CPPUTEST_DIR}/include/CppUTest/MemoryLeakDetectorMallocMacros.h")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -include ${CPPUTEST_DIR}/include/CppUTest/MemoryLeakDetectorNewMacros.h")

find_package(HDF5 COMPONENTS C HL REQUIRED)

set(SRC_LIST
    main.cpp
    tests1.cpp
)

include_directories("${AMICI_DIR}")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories("${MODEL_DIR}")
include_directories("${HDF5_INCLUDE_DIRS}")
include_directories("${AMICI_DIR}/sundials/include")
include_directories("${SUITESPARSE_DIR}/KLU/Include")
include_directories("${SUITESPARSE_DIR}/AMD/Include")
include_directories("${SUITESPARSE_DIR}/SuiteSparse_config")
include_directories("${SUITESPARSE_DIR}/COLAMD/Include")
include_directories("${SUITESPARSE_DIR}/BTF/Include")
include_directories(${CPPUTEST_INCLUDE_DIR})

add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}
${AMICI_DIR}/models/model_dirac/build/libmodel_dirac.a
${CPPUTEST_LIBRARY}
"${SUNDIALS_LIB_DIR}/libsundials_nvecserial.so"
"${SUNDIALS_LIB_DIR}/libsundials_cvodes.so"
"${SUITESPARSE_LIB_DIR}/lib/libcolamd.so"
"${SUITESPARSE_LIB_DIR}/KLU/Lib/libklu.a"
"${SUITESPARSE_LIB_DIR}/BTF/Lib/libbtf.a"
"${SUITESPARSE_LIB_DIR}/AMD/Lib/libamd.a"
"${SUITESPARSE_LIB_DIR}/COLAMD/Lib/libcolamd.a"
"${SUITESPARSE_LIB_DIR}/SuiteSparse_config/libsuitesparseconfig.a"
"${HDF5_HL_LIBRARIES}"
"${HDF5_C_LIBRARIES}"
"-lpthread -ldl -lz -lsz"
"-lm"
)
