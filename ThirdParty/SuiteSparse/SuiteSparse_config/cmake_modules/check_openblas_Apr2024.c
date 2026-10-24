//------------------------------------------------------------------------------
// SuiteSparse/cmake_modules/check_openblas_Apr2024.c
//------------------------------------------------------------------------------

// Copyright (c) 2017-2025, Timothy A. Davis.  All Rights Reserved.
// SPDX-License-Identifier: BSD-3-clause

//------------------------------------------------------------------------------

// check the properties of OpenBLAS 0.3.27 or later (released Apr 4, 2024)

#include <stdio.h>
int openblas_set_num_threads_local (int nthreads) ;
int openblas_get_num_threads (void) ;
int main (void)
{
    int nthreads = openblas_get_num_threads ( ) ;
    int prior = openblas_set_num_threads_local (2) ;
    printf ("%d", nthreads) ;
    return (0) ;
}

