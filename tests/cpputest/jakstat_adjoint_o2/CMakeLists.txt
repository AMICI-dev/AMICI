project(model_jakstat_adjoint_o2_test)

set(MODEL_DIR "${AMICI_DIR}/models/model_jakstat_adjoint_o2")

set(SRC_LIST
    main.cpp
    tests1.cpp
    ../testfunctions.cpp
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories("${MODEL_DIR}")

add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}
    ${AMICI_DIR}/models/model_jakstat_adjoint_o2/build/libmodel_jakstat_adjoint_o2.a
    ${CPPUTEST_LIBRARY}
    amici
    )
    
if($ENV{ENABLE_GCOV_COVERAGE})
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 --coverage")
endif()
    
add_test(NAME jakstat_adjoint_o2 COMMAND ./model_jakstat_adjoint_o2_test -c)
