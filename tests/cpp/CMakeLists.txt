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

# Generate test models
add_custom_target(
  import_test_models
  COMMAND "${CMAKE_SOURCE_DIR}/venv/bin/python" "-c" "from amici.testing.models import import_test_models; import_test_models()"
  VERBATIM
)

# models depend on Upstream::amici
add_library(Upstream::amici ALIAS amici)

# Amici testing library
add_library(amici-testing testfunctions.cpp)
target_compile_definitions(
  amici-testing
  PUBLIC NEW_OPTION_FILE="${CMAKE_CURRENT_SOURCE_DIR}/testOptions.h5"
  PUBLIC HDFFILE="${CMAKE_CURRENT_SOURCE_DIR}/expectedResults.h5"
  PUBLIC HDFFILEWRITE="${CMAKE_CURRENT_SOURCE_DIR}/writeResults.h5")
target_include_directories(amici-testing PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(
  amici-testing
  PUBLIC Upstream::amici
  PUBLIC gtest_main)

# Names of models for which tests are to be run
set(TEST_MODELS
    jakstat_adjoint_o2
    neuron_o2
    )

if(ENABLE_SWIG AND ENABLE_PYTHON)
  add_custom_target(
    python-tests
    COMMAND ${CMAKE_SOURCE_DIR}/scripts/run-python-tests.sh
    DEPENDS)
endif()

add_subdirectory(unittests)

include(ExternalProject)



# Amici testing library
add_library(amici-testing-py testfunctions.cpp)
target_compile_definitions(
  amici-testing-py
  PUBLIC NEW_OPTION_FILE="${CMAKE_CURRENT_SOURCE_DIR}/testOptions.h5"
  PUBLIC HDFFILE="${CMAKE_CURRENT_SOURCE_DIR}/expected_results_py.h5"
  PUBLIC HDFFILEWRITE="${CMAKE_CURRENT_SOURCE_DIR}/writeResults.h5")
target_include_directories(amici-testing-py PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(
  amici-testing-py
  PUBLIC Upstream::amici
  PUBLIC gtest_main)



set(PY_TEST_MODELS
  dirac_py
  events_py
  calvetti_py
  neuron_py
  nested_events_py
  steadystate_py
  robertson_py
  jakstat_adjoint_py
  # TODO jakstat_adjoint_o2 -- requires 2nd-order sensitivities
  # TODO neuron_o2 -- requires 2nd-order sensitivities
)
foreach(MODEL IN ITEMS ${PY_TEST_MODELS})
  # Build model
  string(
    CONCAT MODEL_LIBRARY_DIR
           "${CMAKE_CURRENT_BINARY_DIR}/external_model_${MODEL}-prefix/src/"
           "external_model_${MODEL}-build/")
  string(
    CONCAT MODEL_LIBRARY
           "${MODEL_LIBRARY_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}model_${MODEL}"
           "${CMAKE_STATIC_LIBRARY_SUFFIX}")
  ExternalProject_Add(
    external_model_${MODEL}
    PREFIX ""
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/models/model_${MODEL}/"
    INSTALL_COMMAND ""
    TEST_COMMAND ""
    BUILD_ALWAYS 1
    DEPENDS amici
    CMAKE_ARGS "-DCMAKE_BUILD_TYPE=Debug"
    BUILD_BYPRODUCTS "${MODEL_LIBRARY}")
  # Rebuild if amici files are updated
  ExternalProject_Add_StepDependencies(external_model_${MODEL} build amici)
  # add_dependencies(external_model_${MODEL} import_test_models)
  add_library(model_${MODEL} STATIC IMPORTED)
  add_dependencies(model_${MODEL} external_model_${MODEL})

  set_target_properties(
    model_${MODEL}
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
               "${CMAKE_SOURCE_DIR}/models/model_${MODEL}/" IMPORTED_LOCATION
                                                            "${MODEL_LIBRARY}")
  # Build tests for this model
  # strip "_py" suffix
  string(REPLACE "_py" "" subdir ${MODEL})
  add_subdirectory(${subdir})
  if(TARGET python-tests)
    add_dependencies(python-tests external_model_${MODEL})
  endif()
endforeach()


foreach(MODEL IN ITEMS ${TEST_MODELS})
  # Build model
  string(
    CONCAT MODEL_LIBRARY_DIR
           "${CMAKE_CURRENT_BINARY_DIR}/external_model_${MODEL}-prefix/src/"
           "external_model_${MODEL}-build/")
  string(
    CONCAT MODEL_LIBRARY
           "${MODEL_LIBRARY_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}model_${MODEL}"
           "${CMAKE_STATIC_LIBRARY_SUFFIX}")
  ExternalProject_Add(
    external_model_${MODEL}
    PREFIX ""
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/models/model_${MODEL}/"
    INSTALL_COMMAND ""
    TEST_COMMAND ""
    BUILD_ALWAYS 1
    DEPENDS amici
    CMAKE_ARGS "-DCMAKE_BUILD_TYPE=Debug"
    BUILD_BYPRODUCTS "${MODEL_LIBRARY}")
  # Rebuild if amici files are updated
  ExternalProject_Add_StepDependencies(external_model_${MODEL} build amici)

  add_library(model_${MODEL} STATIC IMPORTED)
  add_dependencies(model_${MODEL} external_model_${MODEL})

  set_target_properties(
    model_${MODEL}
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
               "${CMAKE_SOURCE_DIR}/models/model_${MODEL}/" IMPORTED_LOCATION
                                                            "${MODEL_LIBRARY}")
  # Build tests for this model
  add_subdirectory(${MODEL})
  if(TARGET python-tests)
    add_dependencies(python-tests external_model_${MODEL})
  endif()
endforeach()
