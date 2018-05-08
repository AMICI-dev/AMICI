find_package(PythonInterp 3.6 REQUIRED)

# Configure setup.py and copy to output directory
SET(SETUP_PY_IN ${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in)
SET(SETUP_PY_OUT ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
CONFIGURE_FILE(${SETUP_PY_IN} ${SETUP_PY_OUT})
#CONFIGURE_FILE(paths.py.in ${CMAKE_CURRENT_BINARY_DIR}/amici/paths.py)

# add targets for installation and package generation
add_custom_target(install-python
    COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY_OUT} install --prefix= --user)

add_custom_target(python-wheel
    COMMAND ${PYTHON_EXECUTABLE} setup.py bdist_wheel --dist-dir=${CMAKE_CURRENT_BINARY_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})


# Copy further amici-python-module files to binary_dir
add_custom_command(
        TARGET python-wheel PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/python/amici ${CMAKE_CURRENT_BINARY_DIR}/amici
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR}/amici/src
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/include ${CMAKE_CURRENT_BINARY_DIR}/amici/include
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/swig ${CMAKE_CURRENT_BINARY_DIR}/amici/swig
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/ThirdParty ${CMAKE_CURRENT_BINARY_DIR}/amici/ThirdParty
        COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/README.md ${CMAKE_CURRENT_BINARY_DIR}/
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/amici/libs/
        COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_BINARY_DIR}/libamici.a" ${CMAKE_CURRENT_BINARY_DIR}/amici/libs/
        COMMAND ${CMAKE_COMMAND} -E copy "${PROJECT_SOURCE_DIR}/ThirdParty/sundials/build/lib/libsundials_nvecserial.a" ${CMAKE_CURRENT_BINARY_DIR}/amici/libs/
        COMMAND ${CMAKE_COMMAND} -E copy "${PROJECT_SOURCE_DIR}/ThirdParty/sundials/build/lib/libsundials_cvodes.a" ${CMAKE_CURRENT_BINARY_DIR}/amici/libs/
        COMMAND ${CMAKE_COMMAND} -E copy "${PROJECT_SOURCE_DIR}/ThirdParty/sundials/build/lib/libsundials_idas.a" ${CMAKE_CURRENT_BINARY_DIR}/amici/libs/
        COMMAND ${CMAKE_COMMAND} -E copy "${PROJECT_SOURCE_DIR}/ThirdParty/SuiteSparse/KLU/Lib/libklu.a" ${CMAKE_CURRENT_BINARY_DIR}/amici/libs/
        COMMAND ${CMAKE_COMMAND} -E copy "${PROJECT_SOURCE_DIR}/ThirdParty/SuiteSparse/COLAMD/Lib/libcolamd.a" ${CMAKE_CURRENT_BINARY_DIR}/amici/libs/
        COMMAND ${CMAKE_COMMAND} -E copy "${PROJECT_SOURCE_DIR}/ThirdParty/SuiteSparse/BTF/Lib/libbtf.a" ${CMAKE_CURRENT_BINARY_DIR}/amici/libs/
        COMMAND ${CMAKE_COMMAND} -E copy "${PROJECT_SOURCE_DIR}/ThirdParty/SuiteSparse/AMD/Lib/libamd.a" ${CMAKE_CURRENT_BINARY_DIR}/amici/libs/
        COMMAND ${CMAKE_COMMAND} -E copy "${PROJECT_SOURCE_DIR}/ThirdParty/SuiteSparse/SuiteSparse_config/libsuitesparseconfig.a" ${CMAKE_CURRENT_BINARY_DIR}/amici/libs/
)
#TODO: blas ?

