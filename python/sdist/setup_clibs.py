"""
Provides setuptools clibs for AMICI core, Sundials and SuiteSparse

We could compile all source together into the AMICI base Python module,
however we want to keep the static libs to avoid recompilation for
AMICI-generated models.
"""

import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

PackageInfo = Dict[str, List[Union[str, Tuple[str, Any]]]]
Library = Tuple[str, PackageInfo]


def get_suite_sparse_include_dirs(ss_base_dir: Path) -> List[str]:
    """Suite sparse include directories

    :param ss_base_dir: SuiteSparse base directory
    """
    return list(map(str, [
        ss_base_dir / 'KLU' / 'Include',
        ss_base_dir / 'AMD' / 'Include',
        ss_base_dir / 'COLAMD' / 'Include',
        ss_base_dir / 'BTF' / 'Include',
        ss_base_dir / 'SuiteSparse_config',
        ss_base_dir / 'include',
    ]))


def get_sundials_include_dirs(sundials_base_dir: Path) -> List[str]:
    """Sundials include directories

    :param sundials_base_dir: sundials base directory
    """
    return list(map(str, [
        sundials_base_dir / 'include',
        sundials_base_dir / 'src',
    ]))


def get_sundials_sources(sundials_base_dir: Path) -> List[str]:
    """Get list of Sundials source files

    :param sundials_base_dir: sundials base directory
    """
    return list(map(str, [
        sundials_base_dir / 'src' / 'sunmatrix' / 'dense'
        / 'sunmatrix_dense.c',
        sundials_base_dir / 'src' / 'sunmatrix' / 'band'
        / 'sunmatrix_band.c',
        sundials_base_dir / 'src' / 'sunmatrix' / 'sparse'
        / 'sunmatrix_sparse.c',
        sundials_base_dir / 'src' / 'sunlinsol' / 'spgmr'
        / 'sunlinsol_spgmr.c',
        sundials_base_dir / 'src' / 'sunlinsol' / 'sptfqmr'
        / 'sunlinsol_sptfqmr.c',
        sundials_base_dir / 'src' / 'sunlinsol' / 'klu' / 'sunlinsol_klu.c',
        sundials_base_dir / 'src' / 'sunlinsol' / 'dense'
        / 'sunlinsol_dense.c',
        sundials_base_dir / 'src' / 'sunlinsol' / 'spfgmr'
        / 'sunlinsol_spfgmr.c',
        sundials_base_dir / 'src' / 'sunlinsol' / 'pcg' / 'sunlinsol_pcg.c',
        sundials_base_dir / 'src' / 'sunlinsol' / 'spbcgs'
        / 'sunlinsol_spbcgs.c',
        sundials_base_dir / 'src' / 'sunlinsol' / 'band' / 'sunlinsol_band.c',
        sundials_base_dir / 'src' / 'idas' / 'idas_direct.c',
        sundials_base_dir / 'src' / 'idas' / 'idaa.c',
        sundials_base_dir / 'src' / 'idas' / 'idas_ic.c',
        sundials_base_dir / 'src' / 'idas' / 'idas_nls_stg.c',
        sundials_base_dir / 'src' / 'idas' / 'idas.c',
        sundials_base_dir / 'src' / 'idas' / 'idas_bbdpre.c',
        sundials_base_dir / 'src' / 'idas' / 'idas_spils.c',
        sundials_base_dir / 'src' / 'idas' / 'idas_nls.c',
        sundials_base_dir / 'src' / 'idas' / 'idas_ls.c',
        sundials_base_dir / 'src' / 'idas' / 'idas_io.c',
        sundials_base_dir / 'src' / 'idas' / 'idas_nls_sim.c',
        sundials_base_dir / 'src' / 'idas' / 'idaa_io.c',
        sundials_base_dir / 'src' / 'sundials' / 'sundials_math.c',
        sundials_base_dir / 'src' / 'sundials' / 'sundials_matrix.c',
        sundials_base_dir / 'src' / 'sundials' / 'sundials_direct.c',
        sundials_base_dir / 'src' / 'sundials'
        / 'sundials_nvector_senswrapper.c',
        sundials_base_dir / 'src' / 'sundials' / 'sundials_dense.c',
        sundials_base_dir / 'src' / 'sundials' / 'sundials_nvector.c',
        sundials_base_dir / 'src' / 'sundials' / 'sundials_version.c',
        sundials_base_dir / 'src' / 'sundials' / 'sundials_iterative.c',
        sundials_base_dir / 'src' / 'sundials' / 'sundials_nonlinearsolver.c',
        sundials_base_dir / 'src' / 'sundials' / 'sundials_linearsolver.c',
        sundials_base_dir / 'src' / 'sundials' / 'sundials_band.c',
        sundials_base_dir / 'src' / 'sundials' / 'sundials_futils.c',
        sundials_base_dir / 'src' / 'sunnonlinsol' / 'newton'
        / 'sunnonlinsol_newton.c',
        sundials_base_dir / 'src' / 'sunnonlinsol' / 'fixedpoint'
        / 'sunnonlinsol_fixedpoint.c',
        sundials_base_dir / 'src' / 'nvector' / 'serial' / 'nvector_serial.c',
        sundials_base_dir / 'src' / 'cvodes' / 'cvodes_spils.c',
        sundials_base_dir / 'src' / 'cvodes' / 'cvodes_nls_stg.c',
        sundials_base_dir / 'src' / 'cvodes' / 'cvodes_ls.c',
        sundials_base_dir / 'src' / 'cvodes' / 'cvodes_nls_stg1.c',
        sundials_base_dir / 'src' / 'cvodes' / 'cvodes_bbdpre.c',
        sundials_base_dir / 'src' / 'cvodes' / 'cvodes.c',
        sundials_base_dir / 'src' / 'cvodes' / 'cvodes_bandpre.c',
        sundials_base_dir / 'src' / 'cvodes' / 'cvodea.c',
        sundials_base_dir / 'src' / 'cvodes' / 'cvodes_nls_sim.c',
        sundials_base_dir / 'src' / 'cvodes' / 'cvodea_io.c',
        sundials_base_dir / 'src' / 'cvodes' / 'cvodes_nls.c',
        sundials_base_dir / 'src' / 'cvodes' / 'cvodes_diag.c',
        sundials_base_dir / 'src' / 'cvodes' / 'cvodes_io.c',
        sundials_base_dir / 'src' / 'cvodes' / 'cvodes_direct.c'
    ]))


def get_suite_sparse_sources(ss_base_dir: Path) -> List[str]:
    """Get list of SuiteSparse source files

    :param ss_base_dir: SuiteSparse base directory
    """
    return list(map(str, [
        ss_base_dir / 'KLU' / 'Source' / 'klu_analyze_given.c',
        ss_base_dir / 'KLU' / 'Source' / 'klu_analyze.c',
        ss_base_dir / 'KLU' / 'Source' / 'klu_defaults.c',
        ss_base_dir / 'KLU' / 'Source' / 'klu_diagnostics.c',
        ss_base_dir / 'KLU' / 'Source' / 'klu_dump.c',
        ss_base_dir / 'KLU' / 'Source' / 'klu_extract.c',
        ss_base_dir / 'KLU' / 'Source' / 'klu_factor.c',
        ss_base_dir / 'KLU' / 'Source' / 'klu_free_numeric.c',
        ss_base_dir / 'KLU' / 'Source' / 'klu_free_symbolic.c',
        ss_base_dir / 'KLU' / 'Source' / 'klu_kernel.c',
        ss_base_dir / 'KLU' / 'Source' / 'klu_memory.c',
        ss_base_dir / 'KLU' / 'Source' / 'klu_refactor.c',
        ss_base_dir / 'KLU' / 'Source' / 'klu_scale.c',
        ss_base_dir / 'KLU' / 'Source' / 'klu_sort.c',
        ss_base_dir / 'KLU' / 'Source' / 'klu_solve.c',
        ss_base_dir / 'KLU' / 'Source' / 'klu_tsolve.c',
        ss_base_dir / 'KLU' / 'Source' / 'klu.c',
        ss_base_dir / 'AMD' / 'Source' / 'amd_1.c',
        ss_base_dir / 'AMD' / 'Source' / 'amd_2.c',
        ss_base_dir / 'AMD' / 'Source' / 'amd_aat.c',
        ss_base_dir / 'AMD' / 'Source' / 'amd_control.c',
        ss_base_dir / 'AMD' / 'Source' / 'amd_defaults.c',
        ss_base_dir / 'AMD' / 'Source' / 'amd_dump.c',
        ss_base_dir / 'AMD' / 'Source' / 'amd_global.c',
        ss_base_dir / 'AMD' / 'Source' / 'amd_info.c',
        ss_base_dir / 'AMD' / 'Source' / 'amd_order.c',
        ss_base_dir / 'AMD' / 'Source' / 'amd_post_tree.c',
        ss_base_dir / 'AMD' / 'Source' / 'amd_postorder.c',
        ss_base_dir / 'AMD' / 'Source' / 'amd_preprocess.c',
        ss_base_dir / 'AMD' / 'Source' / 'amd_valid.c',
        ss_base_dir / 'COLAMD' / 'Source' / 'colamd.c',
        ss_base_dir / 'BTF' / 'Source' / 'btf_maxtrans.c',
        ss_base_dir / 'BTF' / 'Source' / 'btf_order.c',
        ss_base_dir / 'BTF' / 'Source' / 'btf_strongcomp.c',
        ss_base_dir / 'SuiteSparse_config' / 'SuiteSparse_config.c',
    ]))


def get_amici_base_sources(
        base_dir: Path,
        with_hdf5: bool = True
) -> List[str]:
    """Get list of source files for the amici base library

    Expects that we are inside $AMICI_ROOT/python/sdist

    Arguments:
        base_dir: AMICI base dir containing ``src/`` and ``include/``
        with_hdf5: compile with HDF5 support
    """
    amici_base_sources = (base_dir / 'src').glob('*.cpp')
    amici_base_sources = [
        str(src) for src in amici_base_sources
        if not re.search(r'(matlab)|(\.(ODE_)?template\.)', str(src))
    ]

    if not with_hdf5:
        hdf5_cpp = base_dir / 'src' / 'hdf5.cpp'
        try:
            # sometimes this fails for unknown reasons...
            amici_base_sources.remove(str(hdf5_cpp))
        except ValueError:
            print(f'Warning: could not find {hdf5_cpp} in '
                  f'{amici_base_sources}')

    return amici_base_sources


def get_lib_sundials(
        sundials_base_dir: Path,
        suitesparse_base_dir: Path,
        extra_compiler_flags: Optional[List[str]] = None
) -> Library:
    """
    Get sundials library build info for setuptools

    :param extra_compiler_flags:
        Extra compiler flags
    """
    if extra_compiler_flags is None:
        extra_compiler_flags = []

    return (
        'sundials',
        {
            'sources': get_sundials_sources(sundials_base_dir),
            'include_dirs': [
                *get_sundials_include_dirs(sundials_base_dir),
                *get_suite_sparse_include_dirs(suitesparse_base_dir),
            ],
            'cflags': [*extra_compiler_flags],
            'cflags_mingw32': ['-Wno-misleading-indentation'],
            'cflags_unix': ['-Wno-misleading-indentation'],
        }
    )


def get_lib_suite_sparse(
        suitesparse_base_dir: Path,
        extra_compiler_flags: Optional[List[str]] = None
) -> Library:
    """
    Get SuiteSparse library build info for setuptools

    :param extra_compiler_flags:
        Extra compiler flags
    """
    if extra_compiler_flags is None:
        extra_compiler_flags = []

    return (
        'suitesparse',
        {
            'sources': get_suite_sparse_sources(suitesparse_base_dir),
            'include_dirs':
                get_suite_sparse_include_dirs(suitesparse_base_dir),
            'cflags': [*extra_compiler_flags],
            'cflags_mingw32': ['-Wno-unused-but-set-variable'],
            'cflags_unix': ['-Wno-unused-but-set-variable']
        })


def get_lib_amici(
        amici_base_dir: Path,
        sundials_base_dir: Path,
        suitesparse_base_dir: Path,
        extra_compiler_flags: List[str] = None,
        h5pkgcfg: Optional[PackageInfo] = None,
        blaspkgcfg: Optional[PackageInfo] = None,
) -> Library:
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
        'sources':
            get_amici_base_sources(
                amici_base_dir,
                with_hdf5=(h5pkgcfg
                           and 'include_dirs' in h5pkgcfg
                           and len(h5pkgcfg['include_dirs']))
            ),
        'include_dirs': [
            str(amici_base_dir / 'include'),
            *get_suite_sparse_include_dirs(suitesparse_base_dir),
            *get_sundials_include_dirs(sundials_base_dir),
            str(amici_base_dir / 'ThirdParty' / 'gsl'),
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
