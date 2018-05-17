find_package(Git)
if(Git_FOUND)
    execute_process(COMMAND sh -c "${GIT_EXECUTABLE} describe --abbrev=4 --dirty=-dirty --always --tags  | tr -d '\n'"
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE AMICI_VERSION
        )
    configure_file(${SRC} ${DST} @ONLY)
endif()

