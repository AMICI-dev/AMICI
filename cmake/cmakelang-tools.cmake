##### Add targets for cmake-format https://cmake-format.readthedocs.io/ ########

# Find all CMakeFiles files
set(all_cxx_source_files
    CMakeLists.txt
    python/CMakeLists.txt
    swig/CMakeLists.txt
    tests/cpp/CMakeLists.txt
    tests/cpp/unittests/CMakeLists.txt
)
list(JOIN all_cxx_source_files " " all_cxx_source_files)

############ cmake-format ############

# Try to find cmake-format and add target if successful
find_program(CMAKE_FORMAT "cmake-format")
if(CMAKE_FORMAT)
    add_custom_target(
        cmake-format
        COMMAND bash -c "${CMAKE_FORMAT} -i ${all_cxx_source_files}"
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        )
else()
    message(STATUS "cmake-format was not found")
endif()
