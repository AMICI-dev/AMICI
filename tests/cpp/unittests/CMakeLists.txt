project(unittests)

find_package(Boost COMPONENTS serialization)

set(SRC_LIST testMisc.cpp testExpData.cpp testSplines.cpp)

add_executable(${PROJECT_NAME} ${SRC_LIST})

target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

if(Boost_SERIALIZATION_FOUND)
  target_sources(${PROJECT_NAME} PRIVATE testSerialization.cpp)
endif()
target_link_libraries(${PROJECT_NAME} amici-testing Upstream::amici
                      Boost::serialization gtest_main)

include(GoogleTest)

gtest_discover_tests(${PROJECT_NAME})
