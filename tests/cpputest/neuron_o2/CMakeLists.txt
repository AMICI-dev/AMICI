project(model_neuron_o2_test)

set(MODEL_DIR "${AMICI_DIR}/models/model_neuron_o2")

set(SRC_LIST
    main.cpp
    tests1.cpp
    ../testfunctions.cpp
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories("${MODEL_DIR}")

add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}
    ${AMICI_DIR}/models/model_neuron_o2/build/libmodel_neuron_o2.a
    ${CPPUTEST_LIBRARY}
    amici
    )
add_test(NAME neuron_o2 COMMAND ./model_neuron_o2_test -c)
