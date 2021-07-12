project(unittests)

find_package(Boost COMPONENTS serialization)

set(SRC_LIST
    tests1.cpp
)

if(Boost_FOUND)
    set(SRC_LIST ${SRC_LIST} testsSerialization.cpp)
    include_directories("${Boost_INCLUDE_DIR}")
endif()

add_executable(${PROJECT_NAME} ${SRC_LIST})

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(${PROJECT_NAME}
    amici-testing
    Upstream::amici
    ${Boost_LIBRARIES}
    gtest_main
    )

include(GoogleTest)

gtest_discover_tests(${PROJECT_NAME})
