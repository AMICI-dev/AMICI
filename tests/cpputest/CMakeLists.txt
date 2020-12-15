project(amiciIntegrationTests)

find_package(CppUTest REQUIRED)
# because Cpputest doesn't seem to care about MEMORY_LEAK_DETECTION=OFF
add_definitions(-DD_MemoryLeakWarningPlugin_h)

# models depend on Upstream::amici
add_library(Upstream::amici ALIAS amici)

set(CMAKE_CXX_FLAGS_OLD "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -include sstream -include functional")

# Amici testing library
add_library(amici-testing testfunctions.cpp)
target_compile_definitions(amici-testing
    PUBLIC NEW_OPTION_FILE="${CMAKE_CURRENT_SOURCE_DIR}/testOptions.h5"
    PUBLIC HDFFILE="${CMAKE_CURRENT_SOURCE_DIR}/expectedResults.h5"
    PUBLIC HDFFILEWRITE="${CMAKE_CURRENT_SOURCE_DIR}/writeResults.h5"
    )
target_include_directories(amici-testing
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
    PUBLIC ${CppUTest_INCLUDE_DIRS}
    )
target_link_libraries(amici-testing
    PUBLIC Upstream::amici
    PUBLIC ${CppUTest_LIBRARIES}
    )

# Names of models for which tests are to be run
set(TEST_MODELS
    dirac
    steadystate
    jakstat_adjoint
    jakstat_adjoint_o2
    neuron neuron_o2
    events
    nested_events
    robertson
    calvetti
    )

if(ENABLE_SWIG AND ENABLE_PYTHON)
    add_custom_target(python-tests
        COMMAND ${CMAKE_SOURCE_DIR}/scripts/run-python-tests.sh
        DEPENDS
    )
endif()

add_subdirectory(unittests)

include(ExternalProject)

foreach(MODEL IN ITEMS ${TEST_MODELS})
    # Build model
    ExternalProject_Add(external_model_${MODEL}
        PREFIX            ""
        SOURCE_DIR        "${CMAKE_SOURCE_DIR}/models/model_${MODEL}/"
        INSTALL_COMMAND   ""
        TEST_COMMAND      ""
        BUILD_ALWAYS      1
        DEPENDS amici
        )
    # Rebuild if amici files are updated
    ExternalProject_Add_StepDependencies(external_model_${MODEL} build amici)

    add_library(model_${MODEL} STATIC IMPORTED)
    add_dependencies(model_${MODEL} external_model_${MODEL})

    set_target_properties(model_${MODEL}
        PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES
        "${CMAKE_SOURCE_DIR}/models/model_${MODEL}/"
        IMPORTED_LOCATION
        "${CMAKE_CURRENT_BINARY_DIR}/external_model_${MODEL}-prefix/src/external_model_${MODEL}-build/${CMAKE_STATIC_LIBRARY_PREFIX}model_${MODEL}${CMAKE_STATIC_LIBRARY_SUFFIX}"
        )
    # Build tests for this model
    add_subdirectory(${MODEL})
    if(TARGET python-tests)
        add_dependencies(python-tests external_model_${MODEL})
    endif()
endforeach()

