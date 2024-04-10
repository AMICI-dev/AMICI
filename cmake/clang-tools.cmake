# Add targets for clang-format and clang-tidy ############

# Find all source files
execute_process(
  COMMAND
    sh -c
    "git ls-tree -r HEAD --name-only src/*.cpp include/amici/*.h | tr '\n' ' '"
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE ALL_CXX_SOURCE_FILES)
# ########### clang-tidy ############

# Try to find clang-format and add target if successful
find_program(CLANG_FORMAT "clang-format")
if(CLANG_FORMAT)
  add_custom_target(
    clang-format
    COMMAND bash -c "${CLANG_FORMAT} -i ${ALL_CXX_SOURCE_FILES}"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
else()
  message(STATUS "clang-format was not found")
endif()

# ########### clang-tidy ############
# Try to find clang-tidy and add target if successful
find_program(CLANG_TIDY "clang-tidy")

if(CLANG_TIDY)
  add_custom_target(
    clang-tidy
    COMMAND
      sh -c
      "${CLANG_TIDY} ${ALL_CXX_SOURCE_FILES} -- -std=c++17 -I${CMAKE_SOURCE_DIR}"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
else()
  message(STATUS "clang-tidy was not found")
endif()
