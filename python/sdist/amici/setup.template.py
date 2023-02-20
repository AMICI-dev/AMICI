"""AMICI model package setup"""
import os

from setuptools import Extension, find_namespace_packages, setup

from amici import amici_path, compiledWithOpenMP, hdf5_enabled, _get_amici_path
from amici.custom_commands import AmiciBuildCMakeExtension
from amici.setuptools import (add_coverage_flags_if_required,
                              add_debug_flags_if_required, add_openmp_flags,
                              get_blas_config, get_hdf5_config)
from cmake_build_extension import CMakeExtension


def get_extension() -> CMakeExtension:
    """Get setuptools extension object for this AMICI model package"""


    # TODO add_coverage_flags_if_required(cxx_flags, linker_flags)
    # TODO add_debug_flags_if_required(cxx_flags, linker_flags)
    # TODO blaspkgcfg = get_blas_config()
    # TODO linker_flags.extend(blaspkgcfg.get('extra_link_args', []))

    # TODO
    # # compiler and linker flags for libamici
    # if 'AMICI_CXXFLAGS' in os.environ:
    #     cxx_flags.extend(os.environ['AMICI_CXXFLAGS'].split(' '))
    # if 'AMICI_LDFLAGS' in os.environ:
    #     linker_flags.extend(os.environ['AMICI_LDFLAGS'].split(' '))

    # ext_include_dirs = [
    #     os.getcwd(),
    #     os.path.join(amici_path, 'include'),
    #     os.path.join(amici_path, "ThirdParty", "gsl"),
    #     os.path.join(amici_path, "ThirdParty", "sundials", "include"),
    #     os.path.join(amici_path, "ThirdParty", "SuiteSparse", "include"),
    #     *h5pkgcfg['include_dirs'],
    #     *blaspkgcfg['include_dirs']
    # ]
    #
    # ext_library_dirs = [
    #     *h5pkgcfg['library_dirs'],
    #     *blaspkgcfg['library_dirs'],
    #     os.path.join(amici_path, 'libs')
    # ]

    # Build shared object
    ext = CMakeExtension(
        name='TPL_MODELNAME._TPL_MODELNAME',
        source_dir='.',
        cmake_configure_options=[
            f"-DAmici_DIR={_get_amici_path()}",
        ],
        # TODO
        # swig_opts=[
        #     '-c++', '-modern', '-outdir', 'TPL_MODELNAME',
        #     '-I%s' % os.path.join(amici_path, 'swig'),
        #     '-I%s' % os.path.join(amici_path, 'include'),
        # ],
        # extra_compile_args=cxx_flags,
        # extra_link_args=linker_flags
    )

    return ext


# Change working directory to setup.py location
os.chdir(os.path.dirname(os.path.abspath(__file__)))

MODEL_EXT = get_extension()

CLASSIFIERS = [
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Science/Research',
    'Operating System :: POSIX :: Linux',
    'Operating System :: MacOS :: MacOS X',
    'Programming Language :: Python',
    'Programming Language :: C++',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
]

CMDCLASS = {
    # CMake-based builds
    'build_ext': AmiciBuildCMakeExtension,
}

# Install
setup(
    name='TPL_MODELNAME',
    cmdclass=CMDCLASS,
    version='TPL_PACKAGE_VERSION',
    description='AMICI-generated module for model TPL_MODELNAME',
    url='https://github.com/AMICI-dev/AMICI',
    author='model-author-todo',
    author_email='model-author-todo',
    ext_modules=[MODEL_EXT],
    packages=find_namespace_packages(),
    install_requires=['amici==TPL_AMICI_VERSION'],
    extras_require={'wurlitzer': ['wurlitzer']},
    python_requires='>=3.8',
    package_data={},
    zip_safe=False,
    classifiers=CLASSIFIERS,
)
