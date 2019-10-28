cmake_minimum_required(VERSION 3.8) # swig_add_library

if(POLICY CMP0078)
  cmake_policy(SET CMP0078 OLD)
endif(POLICY CMP0078)
if(POLICY CMP0074)
  # Use package_ROOT environment variables
  cmake_policy(SET CMP0074 NEW)
endif(POLICY CMP0074)
if(POLICY CMP0086)
  cmake_policy(SET CMP0086 NEW)
endif(POLICY CMP0086)


find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

if(${CMAKE_VERSION} VERSION_LESS "3.12.0")
    find_package(PythonLibs REQUIRED)
    include_directories(${PYTHON_INCLUDE_DIRS})
    set(Python3_LIBRARIES ${PYTHON_LIBRARIES})
else()
    find_package (Python3 COMPONENTS Development)
    include_directories(${Python3_INCLUDE_DIRS})
endif()

set(SWIG_LIBRARY_NAME _${PROJECT_NAME})
set(CMAKE_SWIG_FLAGS "")
set_source_files_properties(${PROJECT_NAME}.i PROPERTIES CPLUSPLUS ON)

# swig does not use INTERFACE_INCLUDE_DIRS of linked libraries, so add manually
get_target_property(AMICI_INCLUDE_DIRS Upstream::amici INTERFACE_INCLUDE_DIRECTORIES)
include_directories(${AMICI_INCLUDE_DIRS} ..)

swig_add_library(${PROJECT_NAME}
    TYPE MODULE
    LANGUAGE python
    SOURCES ${PROJECT_NAME}.i)

swig_link_libraries(${PROJECT_NAME}
    ${Python3_LIBRARIES}
    model)

# configure module setup script
set(SETUP_PY_IN ${Amici_DIR}/model_setup.template.py)
set(SETUP_PY_OUT ${CMAKE_CURRENT_BINARY_DIR}/setup.py)

add_custom_target(install-python
        DEPENDS ${SWIG_LIBRARY_NAME}
        COMMAND python ${SETUP_PY_OUT} install)
