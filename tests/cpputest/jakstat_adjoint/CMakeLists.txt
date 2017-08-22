project(model_jakstat_adjoint_test)

set(MODEL_DIR "${AMICI_DIR}/models/model_jakstat_adjoint")

set(SRC_LIST
    main.cpp
    tests1.cpp
    ../testfunctions.cpp
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories("${MODEL_DIR}")

add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}
    ${AMICI_DIR}/models/model_jakstat_adjoint/build/libmodel_jakstat_adjoint.a
    ${CPPUTEST_LIBRARY}
    amici
    )

add_test(NAME jakstat_adjoint COMMAND ./model_jakstat_adjoint_test -c)
