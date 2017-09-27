project(unittests)

set(SRC_LIST
    main.cpp
    tests1.cpp
    ../testfunctions.cpp
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}
    ${CPPUTEST_LIBRARY}
    amici
    )

if($ENV{ENABLE_GCOV_COVERAGE})
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 --coverage")
endif()

add_test(NAME unittests COMMAND ./unittests -c)
