# ---------------------------------------------------------------
# $Revision: 4957 $
# $Date: 2016-09-23 12:21:47 -0700 (Fri, 23 Sep 2016) $
# ---------------------------------------------------------------
# Programmer:  Steven Smith @ LLNL
# ---------------------------------------------------------------
# LLNS Copyright Start
# Copyright (c) 2014, Lawrence Livermore National Security
# This work was performed under the auspices of the U.S. Department 
# of Energy by Lawrence Livermore National Laboratory in part under 
# Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# LLNS Copyright End
# ---------------------------------------------------------------
# Find KLU library.
# 
# Set library prefixes for Windows
IF(WIN32)
  set(CMAKE_FIND_LIBRARY_PREFIXES "lib")
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib" ".dll")
endif()

### Find include dir
find_path(temp_KLU_INCLUDE_DIR klu.h ${KLU_INCLUDE_DIR})
if (temp_KLU_INCLUDE_DIR)
    set(KLU_INCLUDE_DIR ${temp_KLU_INCLUDE_DIR})
endif()
unset(temp_KLU_INCLUDE_DIR CACHE)
    
if (KLU_LIBRARY)
    # We have (or were given) KLU_LIBRARY - get path to use for other Suitesparse libs
    get_filename_component(KLU_LIBRARY_DIR ${KLU_LIBRARY} PATH)

    # force CACHE update to show user DIR that will be used
    set(KLU_LIBRARY_DIR ${KLU_LIBRARY_DIR} CACHE PATH "" FORCE)
    
else ()
    # find library with user provided directory path
    set(KLU_LIBRARY_NAME klu)
    find_library(KLU_LIBRARY ${KLU_LIBRARY_NAME} ${KLU_LIBRARY_DIR} NO_DEFAULT_PATH)
endif ()
mark_as_advanced(KLU_LIBRARY)

if (NOT AMD_LIBRARY)
    set(AMD_LIBRARY_NAME amd)
    FIND_LIBRARY(AMD_LIBRARY ${AMD_LIBRARY_NAME} ${KLU_LIBRARY_DIR} NO_DEFAULT_PATH)
    mark_as_advanced(AMD_LIBRARY)
endif ()

if (NOT COLAMD_LIBRARY)
    set(COLAMD_LIBRARY_NAME colamd)
    FIND_LIBRARY(COLAMD_LIBRARY ${COLAMD_LIBRARY_NAME} ${KLU_LIBRARY_DIR} NO_DEFAULT_PATH)
    mark_as_advanced(COLAMD_LIBRARY)
endif ()

if (NOT BTF_LIBRARY)
    set(BTF_LIBRARY_NAME btf)
    FIND_LIBRARY( BTF_LIBRARY ${BTF_LIBRARY_NAME} ${KLU_LIBRARY_DIR} NO_DEFAULT_PATH)
    mark_as_advanced(BTF_LIBRARY)
endif ()

if (NOT SUITESPARSECONFIG_LIBRARY)
    set(SUITESPARSECONFIG_LIBRARY_NAME suitesparseconfig)
    FIND_LIBRARY( SUITESPARSECONFIG_LIBRARY ${SUITESPARSECONFIG_LIBRARY_NAME} ${KLU_LIBRARY_DIR} NO_DEFAULT_PATH)
    mark_as_advanced(SUITESPARSECONFIG_LIBRARY)
endif ()

set(KLU_LIBRARIES ${KLU_LIBRARY} ${AMD_LIBRARY} ${COLAMD_LIBRARY} ${BTF_LIBRARY} ${SUITESPARSECONFIG_LIBRARY})
