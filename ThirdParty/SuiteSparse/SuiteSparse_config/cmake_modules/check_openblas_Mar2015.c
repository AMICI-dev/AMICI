//------------------------------------------------------------------------------
// SuiteSparse/cmake_modules/check_openblas_Mar2025.c
//------------------------------------------------------------------------------

// Copyright (c) 2017-2025, Timothy A. Davis.  All Rights Reserved.
// SPDX-License-Identifier: BSD-3-clause

//------------------------------------------------------------------------------

// check the properties of OpenBLAS 0.2.14 or later (released Mar 2015)

#include <stdio.h>
int openblas_get_num_threads (void) ;
int main (void)
{
    int nthreads = openblas_get_num_threads ( ) ;
    printf ("%d", nthreads) ;
    return (0) ;
}

