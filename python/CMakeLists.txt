find_package(PythonInterp 3.6 REQUIRED)

# Configure setup.py and copy to output directory
SET(SETUP_PY_IN ${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in)
SET(SETUP_PY_OUT ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
CONFIGURE_FILE(${SETUP_PY_IN} ${SETUP_PY_OUT})
#CONFIGURE_FILE(paths.py.in ${CMAKE_CURRENT_BINARY_DIR}/amici/paths.py)

# Copy further amici-python-module files to binary_dir
#add_custom_command(
#        TARGET python-wheel PRE_BUILD
execute_process(
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/python/amici ${CMAKE_CURRENT_BINARY_DIR}/amici
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR}/amici/src
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/include ${CMAKE_CURRENT_BINARY_DIR}/amici/include
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/swig ${CMAKE_CURRENT_BINARY_DIR}/amici/swig
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/ThirdParty ${CMAKE_CURRENT_BINARY_DIR}/amici/ThirdParty
        COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/README.md ${CMAKE_CURRENT_BINARY_DIR}/
)

# add targets for installation and package generation
add_custom_target(install-python
    #DEPENDS _amici
    COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY_OUT} install --user)

add_custom_target(python-bdist
    #DEPENDS _amici
    COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY_OUT} bdist)

add_custom_target(python-sdist
    #DEPENDS _amici
    COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY_OUT} sdist)

add_custom_target(python-wheel
    #DEPENDS _amici
    COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY_OUT} bdist_wheel)
