#
# Build additional AMICI interfaces via swig
#

# Use most recent SWIG version available
set(CMAKE_FIND_PACKAGE_SORT_DIRECTION DEC)
find_package(SWIG REQUIRED)

set(SWIG_VERSION_MIN "3.0")
if(${SWIG_VERSION} VERSION_LESS ${SWIG_VERSION_MIN})
  message(FATAL_ERROR "Requiring SWIG>=${SWIG_VERSION_MIN} "
                      "but found only ${SWIG_VERSION}.")
endif()

if(DEFINED ENV{PYTHON_EXECUTABLE})
  set(Python3_EXECUTABLE $ENV{PYTHON_EXECUTABLE})
endif()
if(${CMAKE_VERSION} VERSION_LESS "3.12.0")
  find_package(PythonInterp 3.8 REQUIRED)
else()
  find_package(Python3 COMPONENTS Interpreter Development)
endif()

set(AMICI_INTERFACE_LIST
    ${CMAKE_CURRENT_SOURCE_DIR}/amici.i
    ${CMAKE_CURRENT_SOURCE_DIR}/edata.i
    ${CMAKE_CURRENT_SOURCE_DIR}/rdata.i
    ${CMAKE_CURRENT_SOURCE_DIR}/misc.i
    ${CMAKE_CURRENT_SOURCE_DIR}/model.i
    ${CMAKE_CURRENT_SOURCE_DIR}/model_ode.i
    ${CMAKE_CURRENT_SOURCE_DIR}/model_dae.i
    ${CMAKE_CURRENT_SOURCE_DIR}/solver.i
    ${CMAKE_CURRENT_SOURCE_DIR}/solver_cvodes.i
    ${CMAKE_CURRENT_SOURCE_DIR}/solver_idas.i
    ${CMAKE_CURRENT_SOURCE_DIR}/std_unique_ptr.i
    ${CMAKE_CURRENT_SOURCE_DIR}/hdf5.i
    ${CMAKE_CURRENT_SOURCE_DIR}/abstract_model.i
    ${CMAKE_CURRENT_SOURCE_DIR}/stdvec2numpy.h)

# Add target to show files in IDE
add_custom_target(swigInterface SOURCES ${AMICI_INTERFACE_LIST})

if(AMICI_PYTHON_EXT_ONLY)
    # Build suitesparse and install library and headers
    foreach(subdir SuiteSparse_config BTF AMD CAMD COLAMD KLU)
        execute_process(
          COMMAND make library
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../ThirdParty/SuiteSparse/${subdir}
      )
    endforeach()

    # Build sundials and install library and headers
    execute_process(
        COMMAND ${CMAKE_COMMAND}
        -S${CMAKE_CURRENT_SOURCE_DIR}/../ThirdParty/sundials
        -B${CMAKE_CURRENT_SOURCE_DIR}/../ThirdParty/sundials/build
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_SOURCE_DIR}/../ThirdParty/sundials/build
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} --install ${CMAKE_CURRENT_SOURCE_DIR}/../ThirdParty/sundials/build
    )

    # Get numpy include directory
    # cmake3.19+: COMMAND_ERROR_IS_FATAL
    execute_process(COMMAND ${Python3_EXECUTABLE} -c "from numpy import get_include; print(get_include())"
        OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
        RESULT_VARIABLE RETURN_CODE
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    message("NumPy include directory: ${NUMPY_INCLUDE_DIR}")

    # https://cmake.org/cmake/help/latest/module/UseSWIG.html
    include(${SWIG_USE_FILE})

    # TODO wrap in option
    set(UseSWIG_TARGET_NAME_PREFERENCE STANDARD)
    set_property(SOURCE amici.i PROPERTY CPLUSPLUS ON)
    set_property(SOURCE amici.i PROPERTY OUTPUT_NAME _amici)
    # TODO additional swig flags
    swig_add_library(
      _amici
      LANGUAGE Python
      SOURCES amici.i
      TYPE MODULE)
    set_target_properties(
      _amici
      PROPERTIES
        SWIG_USE_TARGET_INCLUDE_DIRECTORIES TRUE
        INCLUDE_DIRECTORIES
        "${CMAKE_CURRENT_SOURCE_DIR};${CMAKE_CURRENT_SOURCE_DIR}/../include;${NUMPY_INCLUDE_DIR}"
        PREFIX "")

    target_link_libraries(_amici amici Python3::Python)


    # TODO generic extension; follow python naming - _amici.cpython-311-x86_64-linux-gnu.so
    install(
      FILES ${CMAKE_CURRENT_BINARY_DIR}/amici.py
            ${CMAKE_CURRENT_BINARY_DIR}/_amici.so
      DESTINATION .)
endif()
