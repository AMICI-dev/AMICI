#
# Build AMICI python interface
#

# Find Python interpreter, libs and headers
find_package(PythonInterp 3.6 REQUIRED)# TODO to setup.py
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})

set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_BINARY_DIR}/amici)
#set(CMAKE_SWIG_FLAGS "")

set_source_files_properties(
    ${AMICI_SRC_LIST} ${AMICI_DIR}/swig/amici.i
    PROPERTIES CPLUSPLUS ON
    )

SWIG_ADD_LIBRARY(
    amici MODULE
    LANGUAGE python
    SOURCES ${AMICI_DIR}/swig/amici.i
    )

SWIG_LINK_LIBRARIES(amici 
    ${PYTHON_LIBRARIES}
    ${AMICI_DIR}/build/libamici${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUNDIALS_LIB_DIR}/libsundials_nvecserial${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUNDIALS_LIB_DIR}/libsundials_cvodes${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUNDIALS_LIB_DIR}/libsundials_idas${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUITESPARSE_DIR}/KLU/Lib/libklu${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUITESPARSE_DIR}/COLAMD/Lib/libcolamd${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUITESPARSE_DIR}/BTF/Lib/libbtf${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUITESPARSE_DIR}/AMD/Lib/libamd${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SUITESPARSE_DIR}/SuiteSparse_config/libsuitesparseconfig${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${HDF5_HL_LIBRARIES}
    ${HDF5_C_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
    ${BLAS_LIBRARIES}
    -ldl -lz -lm)

# put library into python package
# (SWIG_ADD_LIBRARY directory options don't work in current cmake version)
set_target_properties(_amici
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/amici"
    )

#SET(PYTHON_INSTALL_FILES
#    ${CMAKE_CURRENT_BINARY_DIR}/amici.py
#    $<TARGET_FILE_NAME:_amici>)

# Configure setup.py and copy to output directory
SET(SETUP_PY_IN ${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in)
SET(SETUP_PY_OUT ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
CONFIGURE_FILE(${SETUP_PY_IN} ${SETUP_PY_OUT})
CONFIGURE_FILE(paths.py.in ${CMAKE_CURRENT_BINARY_DIR}/amici/paths.py)

# Copy further amici-python-module files to binary_dir
add_custom_command(
        TARGET _amici POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/python/_amici ${CMAKE_CURRENT_BINARY_DIR}/amici
        COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/README.md ${CMAKE_CURRENT_BINARY_DIR}/
)

# add targets for installation and package generation
add_custom_target(install-python
    DEPENDS _amici
    COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY_OUT} install --user)

add_custom_target(python-bdist
    DEPENDS _amici
    COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY_OUT} bdist)

add_custom_target(python-sdist
    DEPENDS _amici
    COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY_OUT} sdist)

