find_package(PythonInterp 3.6 REQUIRED)

# Configure setup.py and copy to output directory
set(AMICI_VERSION "@AMICI_VERSION@") # to be replaced later
set(SETUP_PY_IN ${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in)
set(SETUP_PY_OUT ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
configure_file(${SETUP_PY_IN} ${SETUP_PY_OUT})

configure_file(${PROJECT_SOURCE_DIR}/python/amici/setup.template.py ${CMAKE_CURRENT_BINARY_DIR}/setup.template.py)

# Copy further amici-python-module files to binary_dir
add_custom_target(python-module-files
    DEPENDS always_rebuild
    COMMENT "Preparing python module directory"
    COMMAND ${CMAKE_COMMAND} -D SRC=${CMAKE_CURRENT_BINARY_DIR}/setup.py
                             -D DST=${CMAKE_CURRENT_BINARY_DIR}/setup.py
                                -P ${PROJECT_SOURCE_DIR}/cmake/configureVersion.cmake
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/amici
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/amici
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/python/amici ${CMAKE_CURRENT_BINARY_DIR}/amici
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR}/amici/src
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/include ${CMAKE_CURRENT_BINARY_DIR}/amici/include
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/swig ${CMAKE_CURRENT_BINARY_DIR}/amici/swig
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/ThirdParty ${CMAKE_CURRENT_BINARY_DIR}/amici/ThirdParty
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/setup.template.py ${CMAKE_CURRENT_BINARY_DIR}/amici/setup.template.py
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

add_custom_target(install-python
    COMMENT "Installing AMICI base python package"
    DEPENDS python-module-files
    COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY_OUT} install --prefix= --user)


# Create python wheel
# Note that we have to run build_ext explicitely before bdist_wheel, otherwise the swig-generated
# amici.py will not be added to the module build folder, because it does not yet exist at the
# build_py stage
add_custom_target(python-wheel
    COMMENT "Creating wheel for AMICI base python package"
    DEPENDS python-module-files
    COMMAND ${PYTHON_EXECUTABLE} setup.py build_ext
    COMMAND ${PYTHON_EXECUTABLE} setup.py bdist_wheel --dist-dir=${CMAKE_CURRENT_BINARY_DIR}
)

add_custom_target(python-sdist
    COMMENT "Creating sdist for AMICI base python package"
    COMMAND ${PYTHON_EXECUTABLE} setup.py sdist --dist-dir=${CMAKE_BINARY_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/sdist
)

add_custom_command(
    OUTPUT always_rebuild
    COMMAND cmake -E echo
    )

