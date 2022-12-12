# --- Add targets for cmake-format https://cmake-format.readthedocs.io/ ---

# Find all CMakeFiles files
set(ALL_CMAKE_FILES
    CMakeLists.txt
    python/CMakeLists.txt
    swig/CMakeLists.txt
    tests/cpp/CMakeLists.txt
    tests/cpp/unittests/CMakeLists.txt
    ${CMAKE_MODULE_PATH}/cmakelang-tools.cmake
    ${CMAKE_MODULE_PATH}/clang-tools.cmake
    ${CMAKE_MODULE_PATH}/version.cmake)
list(JOIN ALL_CMAKE_FILES " " ALL_CMAKE_FILES)

# --- cmake-format ---

# Try to find cmake-format and add target if successful
find_program(CMAKE_FORMAT "cmake-format")
if(CMAKE_FORMAT)
  add_custom_target(
    cmake-format
    COMMAND bash -c "${CMAKE_FORMAT} -i ${ALL_CMAKE_FILES}"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Running cmake-format")
else()
  message(STATUS "cmake-format was not found")
endif()

# --- cmake-lint ---

# Try to find cmake-lint and add target if successful
find_program(CMAKE_LINT "cmake-lint")
if(CMAKE_LINT)
  add_custom_target(
    cmake-lint
    COMMAND bash -c "${CMAKE_LINT} ${ALL_CMAKE_FILES}"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Running cmake-lint")
else()
  message(STATUS "cmake-lint was not found")
endif()
