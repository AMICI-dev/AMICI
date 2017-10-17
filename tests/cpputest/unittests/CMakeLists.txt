project(unittests)

# cannot mix CppuTest new override togeter with Boost
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_OLD}")
find_package(Boost COMPONENTS serialization)

include_directories("${Boost_INCLUDE_DIR}")

set(SRC_LIST
    main.cpp
    tests1.cpp
    ../testfunctions.cpp
)

if(Boost_FOUND)
    set(SRC_LIST ${SRC_LIST} testsSerialization.cpp)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}
    ${CPPUTEST_LIBRARY}
    Upstream::amici
    ${Boost_LIBRARIES}
    )

if($ENV{ENABLE_GCOV_COVERAGE})
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 --coverage")
endif()

add_test(NAME unittests COMMAND ./unittests -c)
