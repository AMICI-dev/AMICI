#-------------------------------------------------------------------------------
# SuiteSparse/SuiteSparse_config/cmake_modules/SuiteSparseLAPACK.cmake
#-------------------------------------------------------------------------------

# SuiteSparse_config, Copyright (c) 2012-2023, Timothy A. Davis.
# All Rights Reserved.
# SPDX-License-Identifier: BSD-3-clause

#-------------------------------------------------------------------------------

# SuiteSparse interface to the Fortran LAPACK library.
# cmake 3.22 is required because BLA_SIZEOF_INTEGER is used.

# The Intel MKL BLAS is highly recommended.  It is free to download (but be
# sure to check their license to make sure you accept it).   See:
# https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.htm

# The use of this package must be preceded with:
# include ( SuiteSparseBLAS )

cmake_minimum_required ( VERSION 3.22 )

if ( BLA_VENDOR STREQUAL "FLAME" )
    # FLAME has the BLAS but not LAPACK

    set ( BLA_VENDOR "Generic" )
    message ( STATUS "Looking for generic LAPACK to use with BLIS/FLAME BLAS" )

    # look for the generic dynamic LAPACK library (usually liblagraph.so)
    find_library ( LAPACK_LIBRARY
        NAMES lapack
        PATH_SUFFIXES lib build )

    # look for the static LAPACK library (usually liblagraph.a)
    if ( MSVC )
        set ( STATIC_SUFFIX .lib )
    else ( )
        set ( STATIC_SUFFIX .a )
    endif ( )
    set ( save ${CMAKE_FIND_LIBRARY_SUFFIXES} )
    set ( CMAKE_FIND_LIBRARY_SUFFIXES ${STATIC_SUFFIX} ${CMAKE_FIND_LIBRARY_SUFFIXES} )
    find_library ( LAPACK_STATIC
        NAMES lapack
        PATH_SUFFIXES lib build)
    set ( CMAKE_FIND_LIBRARY_SUFFIXES ${save} )

    set ( LAPACK_LIBRARIES ${LAPACK_LIBRARY} )

    include (FindPackageHandleStandardArgs)

    find_package_handle_standard_args ( LAPACK
        REQUIRED_VARS LAPACK_LIBRARY
    )

    mark_as_advanced (
        LAPACK_LIBRARY
        LAPACK_STATIC
        LAPACK_LIBRARIES
    )

    set ( BLA_VENDOR "FLAME" )

    if ( LAPACK_FOUND )
        message ( STATUS "LAPACK library: ${LAPACK_LIBRARY}" )
        message ( STATUS "LAPACK static:  ${LAPACK_STATIC}" )
    else ( )
        message ( STATUS "LAPACK not found" )
    endif ( )

else ( )
    # all other cases: BLA_VENDOR works fine for LAPACK
    find_package ( LAPACK REQUIRED )
endif ( )

