get_filename_component(MODEL_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)
project(model_${MODEL_NAME}_test)

set(SRC_LIST
    tests1.cpp
)

add_executable(${PROJECT_NAME} ${SRC_LIST})

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_dependencies(${PROJECT_NAME} external_model_${MODEL_NAME})

target_link_libraries(${PROJECT_NAME}
    amici-testing
    model_${MODEL_NAME}
    gtest_main
)

include(GoogleTest)
gtest_discover_tests(${PROJECT_NAME})
