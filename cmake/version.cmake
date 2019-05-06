find_package(Git)
if(Git_FOUND)
    execute_process(COMMAND sh -c "${GIT_EXECUTABLE} describe --abbrev=4 --dirty=-dirty --always --tags  | cut -c 2- | tr -d '\n'"
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE PROJECT_VERSION
        )
endif()

