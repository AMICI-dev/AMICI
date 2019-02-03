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
        os.path.join('src', 'sunmat_dense', 'fsunmatrix_dense.c'),
        os.path.join('src', 'sunmat_dense', 'sunmatrix_dense.c'),
        os.path.join('src', 'sunlinsol_spgmr', 'fsunlinsol_spgmr.c'),
        os.path.join('src', 'sunlinsol_spgmr', 'sunlinsol_spgmr.c'),
        #os.path.join('src', 'nvec_openmp', 'fnvector_openmp.c'),
        #os.path.join('src', 'nvec_openmp', 'nvector_openmp.c'),
        os.path.join('src', 'sunlinsol_sptfqmr', 'fsunlinsol_sptfqmr.c'),
        os.path.join('src', 'sunlinsol_sptfqmr', 'sunlinsol_sptfqmr.c'),
        os.path.join('src', 'sunlinsol_klu', 'sunlinsol_klu.c'),
        os.path.join('src', 'sunlinsol_klu', 'fsunlinsol_klu.c'),
        #os.path.join('src', 'nvec_petsc', 'nvector_petsc.c'),
        #os.path.join('src', 'sunlinsol_lapackdense',
        #             'fsunlinsol_lapackdense.c'),
        #os.path.join('src', 'sunlinsol_lapackdense',
        #             'sunlinsol_lapackdense.c'),
        #os.path.join('src', 'nvec_pthreads', 'fnvector_pthreads.c'),
        #os.path.join('src', 'nvec_pthreads', 'nvector_pthreads.c'),
        #os.path.join('src', 'sunlinsol_lapackband', 'sunlinsol_lapackband.c'),
        #os.path.join('src', 'sunlinsol_lapackband', 'fsunlinsol_lapackband.c'),
        #os.path.join('src', 'sunlinsol_superlumt', 'sunlinsol_superlumt.c'),
        #os.path.join('src', 'sunlinsol_superlumt', 'fsunlinsol_superlumt.c'),
        # os.path.join('src', 'cvode', 'cvode_spils.c'),
        # os.path.join('src', 'cvode', 'cvode_nls.c'),
        # os.path.join('src', 'cvode', 'cvode.c'),
        # os.path.join('src', 'cvode', 'cvode_bandpre.c'),
        # os.path.join('src', 'cvode', 'cvode_ls.c'),
        # os.path.join('src', 'cvode', 'fcmix', 'fcvewt.c'),
        # os.path.join('src', 'cvode', 'fcmix', 'fcvroot.c'),
        # os.path.join('src', 'cvode', 'fcmix', 'fcvbbd.c'),
        # os.path.join('src', 'cvode', 'fcmix', 'fcvnullnonlinsol.c'),
        # os.path.join('src', 'cvode', 'fcmix', 'fcvode.c'),
        # os.path.join('src', 'cvode', 'fcmix', 'fcvnullmatrix.c'),
        # os.path.join('src', 'cvode', 'fcmix', 'fcvsparse.c'),
        # os.path.join('src', 'cvode', 'fcmix', 'fcvjtimes.c'),
        # os.path.join('src', 'cvode', 'fcmix', 'fcvband.c'),
        # os.path.join('src', 'cvode', 'fcmix', 'fcvnulllinsol.c'),
        # os.path.join('src', 'cvode', 'fcmix', 'fcvbp.c'),
        # os.path.join('src', 'cvode', 'fcmix', 'fcvdense.c'),
        # os.path.join('src', 'cvode', 'fcmix', 'fcvpreco.c'),
        # os.path.join('src', 'cvode', 'cvode_direct.c'),
        # os.path.join('src', 'cvode', 'cvode_diag.c'),
        # os.path.join('src', 'cvode', 'cvode_io.c'),
        # os.path.join('src', 'cvode', 'cvode_bbdpre.c'),
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
        # os.path.join('src', 'arkode', 'arkode_butcher_dirk.c'),
        # os.path.join('src', 'arkode', 'arkode_butcher_erk.c'),
        # os.path.join('src', 'arkode', 'arkode_erkstep_io.c'),
        # os.path.join('src', 'arkode', 'arkode_root.c'),
        # os.path.join('src', 'arkode', 'arkode_mristep_io.c'),
        # os.path.join('src', 'arkode', 'arkode_bbdpre.c'),
        # os.path.join('src', 'arkode', 'arkode_erkstep.c'),
        # os.path.join('src', 'arkode', 'arkode_adapt.c'),
        # os.path.join('src', 'arkode', 'arkode_mristep.c'),
        # os.path.join('src', 'arkode', 'arkode_arkstep_io.c'),
        # os.path.join('src', 'arkode', 'arkode_butcher.c'),
        # os.path.join('src', 'arkode', 'arkode_interp.c'),
        # os.path.join('src', 'arkode', 'fcmix', 'farkband.c'),
        # os.path.join('src', 'arkode', 'fcmix', 'farkbp.c'),
        # os.path.join('src', 'arkode', 'fcmix', 'farkewt.c'),
        # os.path.join('src', 'arkode', 'fcmix', 'farkdensemass.c'),
        # os.path.join('src', 'arkode', 'fcmix', 'farksparsemass.c'),
        # os.path.join('src', 'arkode', 'fcmix', 'farkmtimes.c'),
        # os.path.join('src', 'arkode', 'fcmix', 'farksparse.c'),
        # os.path.join('src', 'arkode', 'fcmix', 'farkdense.c'),
        # os.path.join('src', 'arkode', 'fcmix', 'farkpreco.c'),
        # os.path.join('src', 'arkode', 'fcmix', 'farkbandmass.c'),
        # os.path.join('src', 'arkode', 'fcmix', 'farkjtimes.c'),
        # os.path.join('src', 'arkode', 'fcmix', 'farkroot.c'),
        # os.path.join('src', 'arkode', 'fcmix', 'farkadapt.c'),
        # os.path.join('src', 'arkode', 'fcmix', 'farkode.c'),
        # os.path.join('src', 'arkode', 'fcmix', 'farknullmatrix.c'),
        # os.path.join('src', 'arkode', 'fcmix', 'farkbbd.c'),
        # os.path.join('src', 'arkode', 'fcmix', 'farkmasspreco.c'),
        # os.path.join('src', 'arkode', 'fcmix', 'farknullnonlinsol.c'),
        # os.path.join('src', 'arkode', 'fcmix', 'farkexpstab.c'),
        # os.path.join('src', 'arkode', 'fcmix', 'farknulllinsol.c'),
        # os.path.join('src', 'arkode', 'arkode_arkstep.c'),
        # os.path.join('src', 'arkode', 'arkode.c'),
        # os.path.join('src', 'arkode', 'arkode_ls.c'),
        # os.path.join('src', 'arkode', 'arkode_io.c'),
        # os.path.join('src', 'arkode', 'arkode_arkstep_nls.c'),
        # os.path.join('src', 'arkode', 'arkode_bandpre.c'),
        os.path.join('src', 'sundials', 'sundials_math.c'),
        os.path.join('src', 'sundials', 'sundials_mpi.c'),
        os.path.join('src', 'sundials', 'sundials_sptfqmr.c'),
        os.path.join('src', 'sundials', 'sundials_matrix.c'),
        os.path.join('src', 'sundials', 'sundials_pcg.c'),
        os.path.join('src', 'sundials', 'sundials_direct.c'),
        os.path.join('src', 'sundials', 'sundials_spgmr.c'),
        os.path.join('src', 'sundials', 'sundials_spbcgs.c'),
        os.path.join('src', 'sundials', 'sundials_nvector_senswrapper.c'),
        os.path.join('src', 'sundials', 'sundials_dense.c'),
        os.path.join('src', 'sundials', 'sundials_nvector.c'),
        os.path.join('src', 'sundials', 'sundials_version.c'),
        os.path.join('src', 'sundials', 'sundials_spfgmr.c'),
        os.path.join('src', 'sundials', 'sundials_sparse.c'),
        os.path.join('src', 'sundials', 'sundials_iterative.c'),
        os.path.join('src', 'sundials', 'sundials_nonlinearsolver.c'),
        os.path.join('src', 'sundials', 'sundials_linearsolver.c'),
        os.path.join('src', 'sundials', 'sundials_band.c'),
        os.path.join('src', 'sunlinsol_dense', 'fsunlinsol_dense.c'),
        os.path.join('src', 'sunlinsol_dense', 'sunlinsol_dense.c'),
        os.path.join('src', 'sunmat_band', 'sunmatrix_band.c'),
        os.path.join('src', 'sunmat_band', 'fsunmatrix_band.c'),
        #os.path.join('src', 'nvec_par', 'fnvector_parallel.c'),
        #os.path.join('src', 'nvec_par', 'nvector_parallel.c'),
        os.path.join('src', 'sunlinsol_spfgmr', 'fsunlinsol_spfgmr.c'),
        os.path.join('src', 'sunlinsol_spfgmr', 'sunlinsol_spfgmr.c'),
        os.path.join('src', 'sunnonlinsol', 'newton',
                     'fsunnonlinsol_newton.c'),
        os.path.join('src', 'sunnonlinsol', 'newton', 'sunnonlinsol_newton.c'),
        os.path.join('src', 'sunnonlinsol', 'fixedpoint',
                     'fsunnonlinsol_fixedpoint.c'),
        os.path.join('src', 'sunnonlinsol', 'fixedpoint',
                     'sunnonlinsol_fixedpoint.c'),
        # os.path.join('src', 'ida', 'ida_nls.c'),
        # os.path.join('src', 'ida', 'ida_ic.c'),
        # os.path.join('src', 'ida', 'ida_bbdpre.c'),
        # os.path.join('src', 'ida', 'ida_direct.c'),
        # os.path.join('src', 'ida', 'fcmix', 'fidanullmatrix.c'),
        # os.path.join('src', 'ida', 'fcmix', 'fidapreco.c'),
        # os.path.join('src', 'ida', 'fcmix', 'fidadense.c'),
        # os.path.join('src', 'ida', 'fcmix', 'fidabbd.c'),
        # os.path.join('src', 'ida', 'fcmix', 'fidaewt.c'),
        # os.path.join('src', 'ida', 'fcmix', 'fida.c'),
        # os.path.join('src', 'ida', 'fcmix', 'fidasparse.c'),
        # os.path.join('src', 'ida', 'fcmix', 'fidaroot.c'),
        # os.path.join('src', 'ida', 'fcmix', 'fidaband.c'),
        # os.path.join('src', 'ida', 'fcmix', 'fidanullnonlinsol.c'),
        # os.path.join('src', 'ida', 'fcmix', 'fidajtimes.c'),
        # os.path.join('src', 'ida', 'ida_io.c'),
        # os.path.join('src', 'ida', 'ida.c'),
        # os.path.join('src', 'ida', 'ida_spils.c'),
        # os.path.join('src', 'ida', 'ida_ls.c'),
        os.path.join('src', 'sunmat_sparse', 'sunmatrix_sparse.c'),
        os.path.join('src', 'sunmat_sparse', 'fsunmatrix_sparse.c'),
        #os.path.join('src', 'nvec_openmpdev', 'nvector_openmpdev.c'),
        os.path.join('src', 'nvec_ser', 'nvector_serial.c'),
        os.path.join('src', 'nvec_ser', 'fnvector_serial.c'),
        os.path.join('src', 'sunlinsol_pcg', 'sunlinsol_pcg.c'),
        os.path.join('src', 'sunlinsol_pcg', 'fsunlinsol_pcg.c'),
        #os.path.join('src', 'nvec_parhyp', 'nvector_parhyp.c'),
        os.path.join('src', 'sunlinsol_spbcgs', 'fsunlinsol_spbcgs.c'),
        os.path.join('src', 'sunlinsol_spbcgs', 'sunlinsol_spbcgs.c'),
        os.path.join('src', 'sunlinsol_band', 'fsunlinsol_band.c'),
        os.path.join('src', 'sunlinsol_band', 'sunlinsol_band.c'),
        # os.path.join('src', 'kinsol', 'kinsol.c'),
        # os.path.join('src', 'kinsol', 'kinsol_ls.c'),
        # os.path.join('src', 'kinsol', 'kinsol_spils.c'),
        # os.path.join('src', 'kinsol', 'kinsol_io.c'),
        # os.path.join('src', 'kinsol', 'fcmix', 'fkinbbd.c'),
        # os.path.join('src', 'kinsol', 'fcmix', 'fkinnulllinsol.c'),
        # os.path.join('src', 'kinsol', 'fcmix', 'fkinsparse.c'),
        # os.path.join('src', 'kinsol', 'fcmix', 'fkindense.c'),
        # os.path.join('src', 'kinsol', 'fcmix', 'fkinpreco.c'),
        # os.path.join('src', 'kinsol', 'fcmix', 'fkinsol.c'),
        # os.path.join('src', 'kinsol', 'fcmix', 'fkinband.c'),
        # os.path.join('src', 'kinsol', 'fcmix', 'fkinjtimes.c'),
        # os.path.join('src', 'kinsol', 'fcmix', 'fkinnullmatrix.c'),
        # os.path.join('src', 'kinsol', 'kinsol_bbdpre.c'),
        # os.path.join('src', 'kinsol', 'kinsol_direct.c'),
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


def getLibSundials(extra_compiler_flags=None):
    """Get sundials library build info for setuptools

    Arguments:
        extra_compiler_flags: Extra compiler flags
    """
    if extra_compiler_flags is None:
        extra_compiler_flags = []

    libsundials = ('sundials', {
        'sources': getSundialsSources(),
        'include_dirs': ['amici/ThirdParty/sundials/include',
                         'amici/ThirdParty/sundials/src',
                         'amici/ThirdParty/SuiteSparse/KLU/Include/',
                         'amici/ThirdParty/SuiteSparse/AMD/Include/',
                         'amici/ThirdParty/SuiteSparse/COLAMD/Include/',
                         'amici/ThirdParty/SuiteSparse/BTF/Include/',
                         'amici/ThirdParty/SuiteSparse/SuiteSparse_config',
                         'amici/ThirdParty/SuiteSparse/include'],
        'cflags': ['-Wno-misleading-indentation', *extra_compiler_flags]
    })
    return libsundials


def getLibSuiteSparse(extra_compiler_flags=None):
    """Get SuiteSparse library build info for setuptools

    Arguments:
        extra_compiler_flags: Extra compiler flags
    """
    if extra_compiler_flags is None:
        extra_compiler_flags = []

    libsuitesparse = ('suitesparse', {
        'sources': getSuiteSparseSources(),
        'include_dirs': ['amici/ThirdParty/SuiteSparse/KLU/Include/',
                         'amici/ThirdParty/SuiteSparse/AMD/Include/',
                         'amici/ThirdParty/SuiteSparse/COLAMD/Include/',
                         'amici/ThirdParty/SuiteSparse/BTF/Include/',
                         'amici/ThirdParty/SuiteSparse/SuiteSparse_config',
                         'amici/ThirdParty/SuiteSparse/include'
                         ],
        'cflags': ['-Wno-unused-but-set-variable', *extra_compiler_flags]

    })
    return libsuitesparse


def getLibAmici(extra_compiler_flags=None, h5pkgcfg=None, blaspkgcfg=None):
    """Get AMICI core library build info for setuptools

    Arguments:
        extra_compiler_flags: Extra compiler flags
        h5pkgcfg:  hdf5 package info
        blaspkgcfg: blas package info
    """

    if extra_compiler_flags is None:
        extra_compiler_flags = []

    libamici = ('amici', {
        'sources': getAmiciBaseSources(
            withHDF5=(h5pkgcfg and 'include_dirs' in h5pkgcfg and h5pkgcfg['include_dirs'])
            ),
        'include_dirs': ['amici/include',
                         'amici/ThirdParty/SuiteSparse/KLU/Include/',
                         'amici/ThirdParty/SuiteSparse/AMD/Include/',
                         'amici/ThirdParty/SuiteSparse/COLAMD/Include/',
                         'amici/ThirdParty/SuiteSparse/BTF/Include/',
                         'amici/ThirdParty/SuiteSparse/SuiteSparse_config/',
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

    if blaspkgcfg and 'extra_compile_args' in blaspkgcfg:
        libamici[1]['cflags'].extend(blaspkgcfg['extra_compile_args'])

    return libamici
