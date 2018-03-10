#
# Build AMICI python interface
#

FIND_PACKAGE(PythonLibs REQUIRED)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

SET(CMAKE_SWIG_FLAGS "")

SET_SOURCE_FILES_PROPERTIES(${AMICI_DIR}/swig/edata.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_LIBRARY(edata MODULE LANGUAGE python SOURCES ${AMICI_DIR}/swig/edata.i ${AMICI_DIR}/src/edata.cpp)
SWIG_LINK_LIBRARIES(edata ${PYTHON_LIBRARIES})

SET(PYTHON_INSTALL_FILES
        ${CMAKE_CURRENT_BINARY_DIR}/edata.py
        ${CMAKE_CURRENT_BINARY_DIR}/_edata.so)
        
# Configure setup.py and copy to output directory
SET(SETUP_PY_IN ${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in)
SET(SETUP_PY_OUT ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
CONFIGURE_FILE(${SETUP_PY_IN} ${SETUP_PY_OUT})

add_custom_target(install-python
        DEPENDS _edata
        COMMAND python ${SETUP_PY_OUT} install)
