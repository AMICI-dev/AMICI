#-------------------------------------------------------------------------------
# SuiteSparse/SuiteSparse_config/cmake_modules/SuiteSparseBLAS32.cmake
#-------------------------------------------------------------------------------

# SuiteSparse_config, Copyright (c) 2012-2025, Timothy A. Davis.
# All Rights Reserved.
# SPDX-License-Identifier: BSD-3-clause

#-------------------------------------------------------------------------------

# actions taken when a 32-bit BLAS has been found

if ( SUITESPARSE_REQUIRE_BLAS )
    message ( STATUS "Found ${BLA_VENDOR} 32-bit BLAS" )
endif ( )
add_compile_definitions ( BLAS_${BLA_VENDOR} )
if ( BLA_VENDOR MATCHES "Intel10" )
    add_compile_definitions ( BLAS_Intel10 )
endif ( )
set ( SuiteSparse_BLAS_integer "int32_t" )

# determine the BLAS properties (single/multi-thread and set local threads)
include ( SuiteSparse__blas_threading )

