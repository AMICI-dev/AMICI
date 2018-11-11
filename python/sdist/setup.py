"""Setuptools file for creating AMICI module

This file is based on setuptools alone and does not require CMake. 
All sources are compiled anew.

This file expects to be run from within its directory.

Requires:
- swig3.0
- setuptools
- pkgconfig python+executables
- hdf5 libraries and headers
"""

from setuptools import find_packages, setup, Extension
from setuptools.command.build_ext import build_ext
from setuptools.command.sdist import sdist
from setuptools.command.install_lib import install_lib
from setuptools.command.develop import develop

import os
import sys
import glob
import sysconfig
import subprocess
from shutil import copyfile
import numpy as np # for include directory
import setup_clibs  # Must run from within containing directory

from amici import __version__

from amici.setuptools import (
    getBlasConfig,
    getHdf5Config,
    addCoverageFlagsIfRequired,
    addDebugFlagsIfRequired,
    generateSwigInterfaceFiles,
)

# Extra compiler flags
cxx_flags = []
amici_module_linker_flags = []
define_macros = []

blaspkgcfg = getBlasConfig()
amici_module_linker_flags.extend('-l%s' % l for l in blaspkgcfg['libraries'])

h5pkgcfg = getHdf5Config()

if h5pkgcfg['found']:
    # Manually add linker flags. The libraries passed to Extension will
    # end up in front of the clibs in the linker line and not after, where
    # they are required.
    print("HDF5 library found. Building AMICI with HDF5 support.")
    amici_module_linker_flags.extend(
        ['-l%s' % l for l in ['hdf5_hl_cpp', 'hdf5_hl', 'hdf5_cpp', 'hdf5']])
    extension_sources = [
        'amici/amici_wrap.cxx',  # swig interface
    ]
else:
    print("HDF5 library NOT found. Building AMICI WITHOUT HDF5 support.")
    extension_sources = [
        'amici/amici_wrap_without_hdf5.cxx',  # swig interface
    ]

addCoverageFlagsIfRequired(
    cxx_flags,
    amici_module_linker_flags,
)

addDebugFlagsIfRequired(
    cxx_flags,
    amici_module_linker_flags,
)

# compiler and linker flags for libamici
if 'AMICI_CXXFLAGS' in os.environ:
    cxx_flags.extend(os.environ['AMICI_CXXFLAGS'].split(' '))
if 'AMICI_LDFLAGS' in os.environ:
    amici_module_linker_flags.extend(os.environ['AMICI_LDFLAGS'].split(' '))

libamici = setup_clibs.getLibAmici(
    h5pkgcfg=h5pkgcfg, blaspkgcfg=blaspkgcfg, extra_compiler_flags=cxx_flags)
libsundials = setup_clibs.getLibSundials(extra_compiler_flags=cxx_flags)
libsuitesparse = setup_clibs.getLibSuiteSparse(extra_compiler_flags=cxx_flags)

# Build shared object
amici_module = Extension(
    name='amici._amici',
    sources=extension_sources,
    include_dirs=['amici/include',
                  *libsundials[1]['include_dirs'],
                  *libsuitesparse[1]['include_dirs'],
                  *h5pkgcfg['include_dirs'],
                  *blaspkgcfg['include_dirs'],
                  np.get_include()
                  ],
    # Cannot use here, see above
    # libraries=[
    #    'hdf5_hl_cpp', 'hdf5_hl', 'hdf5_cpp', 'hdf5'
    #],
    define_macros=define_macros,
    library_dirs=[
        *h5pkgcfg['library_dirs'],
        *blaspkgcfg['library_dirs'],
        'amici/libs',  # clib target directory
    ],
    extra_compile_args=['-std=c++11', *cxx_flags],
    extra_link_args=amici_module_linker_flags
)


class my_develop(develop):
    """Custom develop to build clibs"""
    def run(self):

        generateSwigInterfaceFiles()

        self.run_command('build')
        develop.run(self)


class my_install_lib(install_lib):
    """Custom install to allow preserving of debug symbols"""
    def run(self):
        """strip debug symbols

        Returns:

        """
        if 'ENABLE_AMICI_DEBUGGING' in os.environ and os.environ['ENABLE_AMICI_DEBUGGING'] == 'TRUE' and sys.platform == 'darwin':
            search_dir = os.path.join(os.getcwd(),self.build_dir,'amici')
            for file in os.listdir(search_dir):
                if file.endswith('.so'):
                    subprocess.run(['dsymutil',os.path.join(search_dir,file),
                                    '-o',os.path.join(search_dir,file + '.dSYM')])


        # Continue with the actual installation
        install_lib.run(self)


class my_build_ext(build_ext):
    """Custom build_ext to allow keeping otherwise temporary static libs"""

    def run(self):
        """Copy the generated clibs to the extensions folder to be included in the wheel

        Returns:

        """

        if not self.dry_run:  # --dry-run
            libraries = []
            build_clib = ''
            if self.distribution.has_c_libraries():
                # get the previously built static libraries
                build_clib = self.get_finalized_command('build_clib')
                libraries = build_clib.get_library_names() or []
                library_dirs = build_clib.build_clib

            # Module build directory where we want to copy the generated libs
            # to
            if self.inplace == 0:
                build_dir = self.build_lib
            else:
                build_dir = os.getcwd()

            target_dir = os.path.join(build_dir, 'amici', 'libs')
            self.mkpath(target_dir)

            # Copy the generated libs
            for lib in libraries:
                libfilenames = glob.glob(
                    '%s%s*%s.*' % (build_clib.build_clib, os.sep, lib)
                )
                assert len(
                    libfilenames) == 1, "Found unexpected number of files: " % libfilenames

                copyfile(libfilenames[0],
                         os.path.join(target_dir, os.path.basename(libfilenames[0])))

        # Continue with the actual extension building
        build_ext.run(self)


class my_sdist(sdist):
    """Custom sdist to run swig and add the interface files to the source distribution

    Could have relied on letting build_ext run swig. However, that would require any user having swig installed
    during package installation. This way we can postpone that until the package is used to compile generated models.
    """

    def run(self):
        """Setuptools entry-point

        Returns:

        """
        self.runSwig()
        self.saveGitVersion()
        sdist.run(self)

    def runSwig(self):
        """Run swig

        Returns:

        """

        if not self.dry_run:  # --dry-run
            # We create two SWIG interfaces, one with HDF5 support, one without
            generateSwigInterfaceFiles()


    def saveGitVersion(self):
        """Create file with extended version string

        This requires git. We assume that whoever creates the sdist will work inside
        a valid git repository.

        Returns:

        """
        with open("amici/git_version.txt", "w") as f:
            sp = subprocess.run(['git', 'describe',
                                 '--abbrev=4', '--dirty=-dirty',
                                 '--always', '--tags'],
                                 stdout=f)
        assert(sp.returncode == 0)


# Readme as long package description to go on PyPi
# (https://pypi.org/project/amici/)
with open("README.md", "r") as fh:
    long_description = fh.read()

# Remove the "-Wstrict-prototypes" compiler option, which isn't valid for
# C++ to fix warnings.
cfg_vars = sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")


def main():
    # Install
    setup(
        name='amici',
        cmdclass={
            'sdist': my_sdist,
            'build_ext': my_build_ext,
            'install_lib': my_install_lib,
            'develop': my_develop,
        },
        version=__version__,
        description='Advanced multi-language Interface to CVODES and IDAS (%s)',
        long_description=long_description,
        long_description_content_type="text/markdown",
        url='https://github.com/ICB-DCM/AMICI',
        author='Fabian Froehlich, Jan Hasenauer, Daniel Weindl and Paul Stapor',
        author_email='fabian_froehlich@hms.harvard.edu',
        license='BSD',
        libraries=[libamici, libsundials, libsuitesparse],
        ext_modules=[amici_module],
        py_modules=['amici/amici',  # the swig interface
                    'amici/amici_without_hdf5',   # the swig interface
                    ],
        packages=find_packages(),
        package_dir={'amici': 'amici'},
        install_requires=['symengine', 'python-libsbml', 'h5py', 'pandas'],
        python_requires='>=3',
        package_data={
            'amici': ['amici/include/amici/*',
                      'src/*template*',
                      'swig/*',
                      'libs/*',
                      'amici.py',
                      'amici_without_hdf5.py',
                      'setup.py.template',
                      ],
        },
        zip_safe=False,
        include_package_data=True,
        exclude_package_data={
            '': ['README.txt'],
        },
        test_suite="tests",
        classifiers=[
            'Development Status :: 3 - Alpha',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: BSD License',
            'Operating System :: POSIX :: Linux',
            'Operating System :: MacOS :: MacOS X',
            'Programming Language :: Python',
            'Programming Language :: C++',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    )

if __name__ == '__main__':
    main()
