"""
Provides setuptools clibs for AMICI core, Sundials and SuiteSparse

We could compile all source together into the AMICI base Python module,
however we want to keep the static libs to avoid recompilation for
AMICI-generated models.
"""

import os
import glob
import re

from typing import Dict, List, Union, Tuple, Any, Optional

PackageInfo = Dict[str, List[Union[str, Tuple[str, Any]]]]
Library = Tuple[str, PackageInfo]

# suite sparse include directories
suite_sparse_include_dirs = [
    'amici/ThirdParty/SuiteSparse/KLU/Include/',
    'amici/ThirdParty/SuiteSparse/AMD/Include/',
    'amici/ThirdParty/SuiteSparse/COLAMD/Include/',
    'amici/ThirdParty/SuiteSparse/BTF/Include/',
    'amici/ThirdParty/SuiteSparse/SuiteSparse_config',
    'amici/ThirdParty/SuiteSparse/include'
]

# sundials include directories
sundials_include_dirs = [
    'amici/ThirdParty/sundials/include',
    'amici/ThirdParty/sundials/src',
]


def get_sundials_sources() -> List[str]:
    """Get list of Sundials source files"""
    srcs = [
        os.path.join('src', 'sunmatrix', 'dense',   'sunmatrix_dense.c'),
        os.path.join('src', 'sunmatrix', 'band',    'sunmatrix_band.c'),
        os.path.join('src', 'sunmatrix', 'sparse',  'sunmatrix_sparse.c'),
        os.path.join('src', 'sunlinsol', 'spgmr',   'sunlinsol_spgmr.c'),
        os.path.join('src', 'sunlinsol', 'sptfqmr', 'sunlinsol_sptfqmr.c'),
        os.path.join('src', 'sunlinsol', 'klu',     'sunlinsol_klu.c'),
        os.path.join('src', 'sunlinsol', 'dense',   'sunlinsol_dense.c'),
        os.path.join('src', 'sunlinsol', 'spfgmr',  'sunlinsol_spfgmr.c'),
        os.path.join('src', 'sunlinsol', 'pcg',     'sunlinsol_pcg.c'),
        os.path.join('src', 'sunlinsol', 'spbcgs',  'sunlinsol_spbcgs.c'),
        os.path.join('src', 'sunlinsol', 'band',    'sunlinsol_band.c'),
        os.path.join('src', 'idas', 'idas_direct.c'),
        os.path.join('src', 'idas', 'idaa.c'),
        os.path.join('src', 'idas', 'idas_ic.c'),
        os.path.join('src', 'idas', 'idas_nls_stg.c'),
        os.path.join('src', 'idas', 'idas.c'),
        os.path.join('src', 'idas', 'idas_bbdpre.c'),
        os.path.join('src', 'idas', 'idas_spils.c'),
        os.path.join('src', 'idas', 'idas_nls.c'),
        os.path.join('src', 'idas', 'idas_ls.c'),
        os.path.join('src', 'idas', 'idas_io.c'),
        os.path.join('src', 'idas', 'idas_nls_sim.c'),
        os.path.join('src', 'idas', 'idaa_io.c'),
        os.path.join('src', 'sundials', 'sundials_math.c'),
        os.path.join('src', 'sundials', 'sundials_matrix.c'),
        os.path.join('src', 'sundials', 'sundials_direct.c'),
        os.path.join('src', 'sundials', 'sundials_nvector_senswrapper.c'),
        os.path.join('src', 'sundials', 'sundials_dense.c'),
        os.path.join('src', 'sundials', 'sundials_nvector.c'),
        os.path.join('src', 'sundials', 'sundials_version.c'),
        os.path.join('src', 'sundials', 'sundials_iterative.c'),
        os.path.join('src', 'sundials', 'sundials_nonlinearsolver.c'),
        os.path.join('src', 'sundials', 'sundials_linearsolver.c'),
        os.path.join('src', 'sundials', 'sundials_band.c'),
        os.path.join('src', 'sundials', 'sundials_futils.c'),
        os.path.join('src', 'sunnonlinsol', 'newton', 'sunnonlinsol_newton.c'),
        os.path.join('src', 'sunnonlinsol', 'fixedpoint',
                     'sunnonlinsol_fixedpoint.c'),
        os.path.join('src', 'nvector', 'serial', 'nvector_serial.c'),
        os.path.join('src', 'cvodes', 'cvodes_spils.c'),
        os.path.join('src', 'cvodes', 'cvodes_nls_stg.c'),
        os.path.join('src', 'cvodes', 'cvodes_ls.c'),
        os.path.join('src', 'cvodes', 'cvodes_nls_stg1.c'),
        os.path.join('src', 'cvodes', 'cvodes_bbdpre.c'),
        os.path.join('src', 'cvodes', 'cvodes.c'),
        os.path.join('src', 'cvodes', 'cvodes_bandpre.c'),
        os.path.join('src', 'cvodes', 'cvodea.c'),
        os.path.join('src', 'cvodes', 'cvodes_nls_sim.c'),
        os.path.join('src', 'cvodes', 'cvodea_io.c'),
        os.path.join('src', 'cvodes', 'cvodes_nls.c'),
        os.path.join('src', 'cvodes', 'cvodes_diag.c'),
        os.path.join('src', 'cvodes', 'cvodes_io.c'),
        os.path.join('src', 'cvodes', 'cvodes_direct.c')
    ]
    return [os.path.join('amici', 'ThirdParty', 'sundials', src)
            for src in srcs]


def get_suite_sparse_sources() -> List[str]:
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
    return [os.path.join('amici', 'ThirdParty', 'SuiteSparse', src)
            for src in srcs]


def get_amici_base_sources(with_hdf5: bool = True) -> List[str]:
    """Get list of source files for the amici base library

    Expects that we are inside $AMICI_ROOT/python/sdist

    Arguments:
        with_hdf5: compile with HDF5 support
    """

    amici_base_sources = glob.glob(os.path.join('amici', 'src', '*.cpp'))
    amici_base_sources = [src for src in amici_base_sources
                          if not re.search(r'(matlab)|(\.(ODE_)?template\.)', src)]

    if not with_hdf5:
        hdf5_cpp = os.path.join('amici', 'src', 'hdf5.cpp')
        try:
            # sometimes this fails for unknwon reasons...
            amici_base_sources.remove(hdf5_cpp)
        except ValueError:
            print(f'Warning: could not find {hdf5_cpp} in '
                  f'{amici_base_sources}')

    return amici_base_sources


def get_lib_sundials(extra_compiler_flags: Optional[List[str]] = None) -> \
        Library:
    """
    Get sundials library build info for setuptools

    :param extra_compiler_flags:
        Extra compiler flags
    """
    if extra_compiler_flags is None:
        extra_compiler_flags = []

    libsundials = ('sundials', {
        'sources': get_sundials_sources(),
        'include_dirs': [*sundials_include_dirs,
                         *suite_sparse_include_dirs,
                         ],
        'cflags': [*extra_compiler_flags],
        'cflags_mingw32': ['-Wno-misleading-indentation'],
        'cflags_unix': ['-Wno-misleading-indentation'],
    })
    return libsundials


def get_lib_suite_sparse(extra_compiler_flags: Optional[List[str]] = None) -> \
        Library:
    """
    Get SuiteSparse library build info for setuptools

    :param extra_compiler_flags:
        Extra compiler flags
    """
    if extra_compiler_flags is None:
        extra_compiler_flags = []

    libsuitesparse = ('suitesparse', {
        'sources': get_suite_sparse_sources(),
        'include_dirs': suite_sparse_include_dirs,
        'cflags': [*extra_compiler_flags],
        'cflags_mingw32': ['-Wno-unused-but-set-variable'],
        'cflags_unix': ['-Wno-unused-but-set-variable']
    })
    return libsuitesparse


def get_lib_amici(extra_compiler_flags: List[str] = None,
                  h5pkgcfg: Optional[PackageInfo] = None,
                  blaspkgcfg: Optional[PackageInfo] = None) -> Library:
    """
    Get AMICI core library build info for setuptools

    :param extra_compiler_flags:
        Extra compiler flags

    :param h5pkgcfg:
        hdf5 package info

    :param blaspkgcfg:
        blas package info

    """

    if extra_compiler_flags is None:
        extra_compiler_flags = []

    libamici = ('amici', {
        'sources': get_amici_base_sources(
            with_hdf5=(h5pkgcfg
                       and 'include_dirs' in h5pkgcfg
                       and len(h5pkgcfg['include_dirs']))
            ),
        'include_dirs': ['amici/include',
                         *suite_sparse_include_dirs,
                         *sundials_include_dirs,
                         'amici/ThirdParty/gsl/',
                         ],
        'cflags': [*extra_compiler_flags],
        'cflags_mingw32': ['-std=c++14'],
        'cflags_unix': ['-std=c++14'],
        'cflags_msvc': ['/std:c++14'],
        'macros': [],
    })

    if h5pkgcfg and 'include_dirs' in h5pkgcfg:
        libamici[1]['include_dirs'].extend(h5pkgcfg['include_dirs'])

    if h5pkgcfg and 'define_macros' in h5pkgcfg:
        libamici[1]['macros'].extend(h5pkgcfg['define_macros'])

    if blaspkgcfg and 'include_dirs' in blaspkgcfg:
        libamici[1]['include_dirs'].extend(blaspkgcfg['include_dirs'])

    if blaspkgcfg and 'define_macros' in blaspkgcfg:
        libamici[1]['macros'].extend(blaspkgcfg['define_macros'])

    if blaspkgcfg and 'extra_compile_args' in blaspkgcfg:
        libamici[1]['cflags'].extend(blaspkgcfg['extra_compile_args'])

    return libamici
