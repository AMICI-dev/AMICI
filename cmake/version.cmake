find_package(Git)
if(Git_FOUND)
  execute_process(
    COMMAND
      sh -c
      "'${GIT_EXECUTABLE}' describe --abbrev=4 --dirty=-dirty --always --tags  | cut -c 2- | tr -d '\n' | sed s/-/./"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE PROJECT_VERSION_GIT)
endif()

# get project root directory
get_filename_component(CMAKE_PARENT_LIST_DIR ${CMAKE_CURRENT_LIST_DIR}
                       DIRECTORY)
file(STRINGS "${CMAKE_PARENT_LIST_DIR}/version.txt" PROJECT_VERSION_RAW)

# try to extract major.minor.patch
# (ignore any Python .dev, .post, etc. suffixes)
string(REGEX MATCH "([0-9]+\\.[0-9]+\\.[0-9]+)" PROJECT_VERSION "${PROJECT_VERSION_RAW}")
if(NOT PROJECT_VERSION)
  message(FATAL_ERROR "Could not extract version number from version.txt")
endif()

message(DEBUG "Version number: ${PROJECT_VERSION}")
