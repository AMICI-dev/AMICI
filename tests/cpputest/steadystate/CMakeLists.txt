get_filename_component(MODEL_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)
project(model_${MODEL_NAME}_test)

set(SRC_LIST
    ../main.cpp
    tests1.cpp
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CppUTest_INCLUDE_DIRS})

add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}
    amici-testing
    model_${MODEL_NAME}
)

add_test(NAME ${PROJECT_NAME} COMMAND ./${PROJECT_NAME} -c)
