//------------------------------------------------------------------------------
// SuiteSparse/cmake_modules/check_mkl.c
//------------------------------------------------------------------------------

// Copyright (c) 2017-2025, Timothy A. Davis.  All Rights Reserved.
// SPDX-License-Identifier: BSD-3-clause

//------------------------------------------------------------------------------

// check if the Intel MKL is single-threaded or multi-threaded

#include <stdio.h>
#include <mkl.h>
int main (void)
{
    int nthreads = mkl_get_max_threads ( ) ;
    int prior = mkl_set_num_threads_local (2) ;
    printf ("%d", nthreads) ;
    return (0) ;
}

