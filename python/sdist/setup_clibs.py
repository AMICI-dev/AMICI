"""Provides setuptools clibs for AMICI core, Sundials and SuiteSparse

We could compile all source together into the AMICI base Python module,
however we want to keep the static libs to avoid recompilation for AMICI-generated models.
"""

import os
import glob
import re


def getSundialsSources():
    """Get list of Sundials source files"""
    srcs = [
        os.path.join('src', 'cvodes', 'cvodes_band.c'),
        os.path.join('src', 'cvodes', 'cvodes_bandpre.c'),
        os.path.join('src', 'cvodes', 'cvodes_bbdpre.c'),
        os.path.join('src', 'cvodes', 'cvodes_direct.c'),
        os.path.join('src', 'cvodes', 'cvodes_dense.c'),
        os.path.join('src', 'cvodes', 'cvodes_sparse.c'),
        os.path.join('src', 'cvodes', 'cvodes_diag.c'),
        os.path.join('src', 'cvodes', 'cvodea.c'),
        os.path.join('src', 'cvodes', 'cvodes.c'),
        os.path.join('src', 'cvodes', 'cvodes_io.c'),
        os.path.join('src', 'cvodes', 'cvodea_io.c'),
        os.path.join('src', 'cvodes', 'cvodes_spils.c'),
        os.path.join('src', 'cvodes', 'cvodes_spbcgs.c'),
        os.path.join('src', 'cvodes', 'cvodes_spgmr.c'),
        os.path.join('src', 'cvodes', 'cvodes_sptfqmr.c'),
        os.path.join('src', 'cvodes', 'cvodes_klu.c'),
        os.path.join('src', 'idas', 'idas.c'),
        os.path.join('src', 'idas', 'idas_sptfqmr.c'),
        os.path.join('src', 'idas', 'idas_spils.c'),
        os.path.join('src', 'idas', 'idas_spgmr.c'),
        os.path.join('src', 'idas', 'idas_spbcgs.c'),
        os.path.join('src', 'idas', 'idas_sparse.c'),
        os.path.join('src', 'idas', 'idas_klu.c'),
        os.path.join('src', 'idas', 'idas_io.c'),
        os.path.join('src', 'idas', 'idas_ic.c'),
        os.path.join('src', 'idas', 'idas_direct.c'),
        os.path.join('src', 'idas', 'idas_dense.c'),
        os.path.join('src', 'idas', 'idas_bbdpre.c'),
        os.path.join('src', 'idas', 'idas_band.c'),
        os.path.join('src', 'idas', 'idaa.c'),
        os.path.join('src', 'idas', 'idaa_io.c'),
        os.path.join('src', 'sundials', 'sundials_band.c'),
        os.path.join('src', 'sundials', 'sundials_dense.c'),
        os.path.join('src', 'sundials', 'sundials_sparse.c'),
        os.path.join('src', 'sundials', 'sundials_iterative.c'),
        os.path.join('src', 'sundials', 'sundials_nvector.c'),
        os.path.join('src', 'sundials', 'sundials_direct.c'),
        os.path.join('src', 'sundials', 'sundials_spbcgs.c'),
        os.path.join('src', 'sundials', 'sundials_spgmr.c'),
        os.path.join('src', 'sundials', 'sundials_sptfqmr.c'),
        os.path.join('src', 'sundials', 'sundials_math.c'),
        os.path.join('src', 'nvec_ser', 'nvector_serial.c'),
    ]
    return [os.path.join('amici', 'ThirdParty', 'sundials', src) for src in srcs]


def getSuiteSparseSources():
    """Get list of SuiteSparse source files"""
    srcs = [
        os.path.join('KLU', 'Source', 'klu_analyze_given.c'),
        os.path.join('KLU', 'Source', 'klu_analyze.c'),
        os.path.join('KLU', 'Source', 'klu_defaults.c'),
        os.path.join('KLU', 'Source', 'klu_diagnostics.c'),
        os.path.join('KLU', 'Source', 'klu_dump.c'),
        os.path.join('KLU', 'Source', 'klu_extract.c'),
        os.path.join('KLU', 'Source', 'klu_factor.c'),
        os.path.join('KLU', 'Source', 'klu_free_numeric.c'),
        os.path.join('KLU', 'Source', 'klu_free_symbolic.c'),
        os.path.join('KLU', 'Source', 'klu_kernel.c'),
        os.path.join('KLU', 'Source', 'klu_memory.c'),
        os.path.join('KLU', 'Source', 'klu_refactor.c'),
        os.path.join('KLU', 'Source', 'klu_scale.c'),
        os.path.join('KLU', 'Source', 'klu_sort.c'),
        os.path.join('KLU', 'Source', 'klu_solve.c'),
        os.path.join('KLU', 'Source', 'klu_tsolve.c'),
        os.path.join('KLU', 'Source', 'klu.c'),
        os.path.join('AMD', 'Source', 'amd_1.c'),
        os.path.join('AMD', 'Source', 'amd_2.c'),
        os.path.join('AMD', 'Source', 'amd_aat.c'),
        os.path.join('AMD', 'Source', 'amd_control.c'),
        os.path.join('AMD', 'Source', 'amd_defaults.c'),
        os.path.join('AMD', 'Source', 'amd_dump.c'),
        os.path.join('AMD', 'Source', 'amd_global.c'),
        os.path.join('AMD', 'Source', 'amd_info.c'),
        os.path.join('AMD', 'Source', 'amd_order.c'),
        os.path.join('AMD', 'Source', 'amd_post_tree.c'),
        os.path.join('AMD', 'Source', 'amd_postorder.c'),
        os.path.join('AMD', 'Source', 'amd_preprocess.c'),
        os.path.join('AMD', 'Source', 'amd_valid.c'),
        os.path.join('COLAMD', 'Source', 'colamd.c'),
        os.path.join('BTF', 'Source', 'btf_maxtrans.c'),
        os.path.join('BTF', 'Source', 'btf_order.c'),
        os.path.join('BTF', 'Source', 'btf_strongcomp.c'),
        os.path.join('SuiteSparse_config', 'SuiteSparse_config.c'),
    ]
    return [os.path.join('amici', 'ThirdParty', 'SuiteSparse', src) for src in srcs]


def getAmiciBaseSources(withHDF5=True):
    """Get list of source files for the amici base library

    Expects that we are inside $AMICI_ROOT/python/sdist
    
    Arguments:
        withHDF5: compile with HDF5 support
    """
    
    amiciBaseSources = glob.glob('amici{s}src{s}*.cpp'.format(s=os.sep))
    amiciBaseSources = [src for src in amiciBaseSources if not re.search(
        r'(matlab)|(\.template\.)', src)]
    
    if not withHDF5:
        try:
            # sometimes this fails for unknwon reasons...
            amiciBaseSources.remove('amici{s}src{s}hdf5.cpp'.format(s=os.sep))
        except ValueError:
            print('Warning: could not find %s in %s' % ('amici{s}src{s}hdf5.cpp'.format(s=os.sep), 
                                                        amiciBaseSources)) 
    
    return amiciBaseSources


def getLibSundials(extra_compiler_flags=[]):
    """Get sundials library build info for setuptools

    Arguments:
        extra_compiler_flags: Extra compiler flags
    """

    libsundials = ('sundials', {
        'sources': getSundialsSources(),
        'include_dirs': ['amici/ThirdParty/sundials/include',
                         'amici/ThirdParty/sundials/src',
                         'amici/ThirdParty/SuiteSparse/KLU/Include/',
                         'amici/ThirdParty/SuiteSparse/AMD/Include/',
                         'amici/ThirdParty/SuiteSparse/COLAMD/Include/',
                         'amici/ThirdParty/SuiteSparse/BTF/Include/',
                         'amici/ThirdParty/SuiteSparse/SuiteSparse_config/Include/',
                         'amici/ThirdParty/SuiteSparse/include'],
        'cflags': ['-Wno-misleading-indentation', *extra_compiler_flags]
    })
    return libsundials


def getLibSuiteSparse(extra_compiler_flags=[]):
    """Get SuiteSparse library build info for setuptools

    Arguments:
        extra_compiler_flags: Extra compiler flags
    """

    libsuitesparse = ('suitesparse', {
        'sources': getSuiteSparseSources(),
        'include_dirs': ['amici/ThirdParty/SuiteSparse/KLU/Include/',
                         'amici/ThirdParty/SuiteSparse/AMD/Include/',
                         'amici/ThirdParty/SuiteSparse/COLAMD/Include/',
                         'amici/ThirdParty/SuiteSparse/BTF/Include/',
                         'amici/ThirdParty/SuiteSparse/SuiteSparse_config/Include/',
                         'amici/ThirdParty/SuiteSparse/include'
                         ],
        'cflags': ['-Wno-unused-but-set-variable', *extra_compiler_flags]

    })
    return libsuitesparse


def getLibAmici(extra_compiler_flags=[], h5pkgcfg=None, blaspkgcfg=None):
    """Get AMICI core library build info for setuptools

    Arguments:
        extra_compiler_flags: Extra compiler flags
        h5pkgcfg:  hdf5 package info
        blaspkgcfg: blas package info
    """

    libamici = ('amici', {
        'sources': getAmiciBaseSources(
            withHDF5=(h5pkgcfg and 'include_dirs' in h5pkgcfg and h5pkgcfg['include_dirs'])
            ),
        'include_dirs': ['amici/include',
                         'amici/ThirdParty/SuiteSparse/KLU/Include/',
                         'amici/ThirdParty/SuiteSparse/AMD/Include/',
                         'amici/ThirdParty/SuiteSparse/COLAMD/Include/',
                         'amici/ThirdParty/SuiteSparse/BTF/Include/',
                         'amici/ThirdParty/SuiteSparse/SuiteSparse_config/Include/',
                         'amici/ThirdParty/SuiteSparse/include',
                         'amici/ThirdParty/sundials/include',
                         'amici/ThirdParty/sundials/src'
                         ],
        'cflags': ['-std=c++11', *extra_compiler_flags]
    })

    if h5pkgcfg and 'include_dirs' in h5pkgcfg:
        libamici[1]['include_dirs'].extend(h5pkgcfg['include_dirs'])

    if blaspkgcfg and 'include_dirs' in blaspkgcfg:
        libamici[1]['include_dirs'].extend(blaspkgcfg['include_dirs'])

    return libamici
