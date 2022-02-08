# ------------------------------------------------------------------------------
# Set up google test
# ------------------------------------------------------------------------------

# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(
    COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download)
if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(
    COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download)
if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker settings on Windows
set(gtest_force_shared_crt
    ON
    CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines the gtest and gtest_main
# targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build EXCLUDE_FROM_ALL)

# ------------------------------------------------------------------------------
# AMICI C++ tests
# ------------------------------------------------------------------------------

project(amiciIntegrationTests)

# models depend on Upstream::amici
add_library(Upstream::amici ALIAS amici)

# Amici testing library
add_library(amici-testing testfunctions.cpp)
target_compile_definitions(amici-testing
    PUBLIC NEW_OPTION_FILE="${CMAKE_CURRENT_SOURCE_DIR}/testOptions.h5"
    PUBLIC HDFFILE="${CMAKE_CURRENT_SOURCE_DIR}/expectedResults.h5"
    PUBLIC HDFFILEWRITE="${CMAKE_CURRENT_SOURCE_DIR}/writeResults.h5"
    )
target_include_directories(amici-testing
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
    )
target_link_libraries(amici-testing
    PUBLIC Upstream::amici
    PUBLIC gtest_main
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
    string(CONCAT MODEL_LIBRARY_DIR
        "${CMAKE_CURRENT_BINARY_DIR}/external_model_${MODEL}-prefix/src/"
        "external_model_${MODEL}-build/")
    string(CONCAT MODEL_LIBRARY
        "${MODEL_LIBRARY_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}model_${MODEL}"
        "${CMAKE_STATIC_LIBRARY_SUFFIX}")
    ExternalProject_Add(external_model_${MODEL}
        PREFIX            ""
        SOURCE_DIR        "${CMAKE_SOURCE_DIR}/models/model_${MODEL}/"
        INSTALL_COMMAND   ""
        TEST_COMMAND      ""
        BUILD_ALWAYS      1
        DEPENDS amici
        BUILD_BYPRODUCTS  "${MODEL_LIBRARY}"
        )
    # Rebuild if amici files are updated
    ExternalProject_Add_StepDependencies(external_model_${MODEL} build amici)

    add_library(model_${MODEL} STATIC IMPORTED)
    add_dependencies(model_${MODEL} external_model_${MODEL})

    set_target_properties(model_${MODEL}
        PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES
        "${CMAKE_SOURCE_DIR}/models/model_${MODEL}/"
        IMPORTED_LOCATION "${MODEL_LIBRARY}"
        )
    # Build tests for this model
    add_subdirectory(${MODEL})
    if(TARGET python-tests)
        add_dependencies(python-tests external_model_${MODEL})
    endif()
endforeach()

