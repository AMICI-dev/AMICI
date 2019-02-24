project(unittests)

# cannot mix CppuTest new override togeter with Boost
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_OLD}")
find_package(Boost COMPONENTS serialization)

set(SRC_LIST
    ../main.cpp
    tests1.cpp
)

if(Boost_FOUND)
    set(SRC_LIST ${SRC_LIST} testsSerialization.cpp)
    include_directories("${Boost_INCLUDE_DIR}")
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}
    amici-testing
    Upstream::amici
    ${Boost_LIBRARIES}
    )

add_test(NAME unittests COMMAND ./unittests -c)
