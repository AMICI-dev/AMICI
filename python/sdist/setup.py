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
import os
import sys
import glob
import sysconfig
from shutil import copyfile

import setup_clibs # Must run from within containing directory

# Extra compiler flags
cxx_flags = []
# TODO: move flexible blas
amici_module_linker_flags = ['-lcblas']

# Find HDF5 include dir and libs
import pkgconfig
h5pkgcfg = pkgconfig.parse("hdf5")

# Manually add linker flags. The libraries passed to Extension will
# end up in front of the clibs in the linker line and not after, where they are required.
amici_module_linker_flags.extend(
    ['-l%s' % l for l in ['hdf5_hl_cpp', 'hdf5_hl', 'hdf5_cpp', 'hdf5']])
if 'ENABLE_GCOV_COVERAGE' in os.environ and os.environ['ENABLE_GCOV_COVERAGE'] == 'TRUE':
    cxx_flags.extend(['-g', '-O0',  '--coverage'])
    amici_module_linker_flags.append('--coverage')

libamici = setup_clibs.getLibAmici(h5pkgcfg=h5pkgcfg, extra_compiler_flags=cxx_flags)
libsundials = setup_clibs.getLibSundials(extra_compiler_flags=cxx_flags)
libsuitesparse = setup_clibs.getLibSuiteSparse(extra_compiler_flags=cxx_flags)


# Build shared object
amici_module = Extension(
    'amici/_amici',
    sources=[
        'amici/amici_wrap.cxx', # swig interface
    ],
    
    include_dirs=['amici/include', 
                  *libsundials[1]['include_dirs'], 
                  *libsuitesparse[1]['include_dirs'],
    #              'amici/ThirdParty/SuiteSparse/KLU/Include/',
    #           'amici/ThirdParty/SuiteSparse/AMD/Include/',
    #              'amici/ThirdParty/SuiteSparse/COLAMD/Include/',
    #             'amici/ThirdParty/SuiteSparse/BTF/Include/',
    #              'amici/ThirdParty/SuiteSparse/SuiteSparse_config/Include/',
    #              'amici/ThirdParty/SuiteSparse/include',
    #              'amici/ThirdParty/sundials/include',
    #             'amici/ThirdParty/sundials/src',
                  *h5pkgcfg['include_dirs']],  # NOTE: requires that pkgconfig knows about hdf5
    # Cannot use here, see above
    #libraries=[
    #    'hdf5_hl_cpp', 'hdf5_hl', 'hdf5_cpp', 'hdf5'
    #],
    library_dirs=[
        *h5pkgcfg['library_dirs'],
        'amici/libs/',
    ],
    extra_compile_args=['-std=c++11', *cxx_flags],
    extra_link_args=amici_module_linker_flags
)


class my_build_ext(build_ext):
    """Custom build_ext to allow keeping otherwise temporary static libs"""
    
    def run(self):
        """Copy the generated clibs to the extensions folder to be included in the wheel"""
        
        if not self.dry_run: # --dry-run
            if self.distribution.has_c_libraries():
                # get the previously built static libraries
                build_clib = self.get_finalized_command('build_clib')
                libraries = build_clib.get_library_names() or []
                library_dirs = build_clib.build_clib

            # Module build directory where we want to copy the generated libs to
            target_dir = os.path.join(self.build_lib, 'amici/libs')
            self.mkpath(target_dir)

            # Copy the generated libs
            for lib in libraries:
                libfilenames = glob.glob('%s/*%s.*' %
                                         (build_clib.build_clib, lib))
                assert(len(libfilenames) == 1)
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
        """Run swig"""
      
        if not self.dry_run: # --dry-run
            import subprocess
            sp = subprocess.run(['swig3.0',
                            '-c++',
                            '-python',
                            '-Iamici/swig', '-Iamici/include', 
                            '-outdir', '%s/amici' % os.path.abspath(os.getcwd()), 
                            '-o', 'amici/amici_wrap.cxx', 
                            'amici/swig/amici.i'])
            assert(sp.returncode == 0)
        sdist.run(self)

# Readme as long package description to go on PyPi (https://pypi.org/project/amici/) 
with open("README.md", "r") as fh:
    long_description = fh.read()

def getVersionNumber():
    return '0.6a4'

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
            'build_ext': my_build_ext
            },
        version=getVersionNumber(),
        description='Advanced multi-language Interface to CVODES and IDAS',
        long_description=long_description,
        long_description_content_type="text/markdown",
        url='https://github.com/ICB-DCM/AMICI',
        author='Fabian Froehlich, Jan Hasenauer, Daniel Weindl and Paul Stapor',
        author_email='fabian.froehlich@helmholtz-muenchen.de',
        license='BSD',
        libraries=[libamici, libsundials, libsuitesparse],
        ext_modules=[amici_module],
        py_modules=['amici/amici'], # the swig interface
        packages=find_packages(),
        package_dir={'amici': 'amici'},
        install_requires=['symengine', 'python-libsbml', 'h5py', 'pkgconfig'],
        python_requires='>=3',
        package_data={
            # TODO Could possibly clean up here, check with MANIFEST.in
            'amici': ['amici/include/amici/*',
                      'src/*template*',
                      'swig/*',
                      'libs/*',
                      'amici.py',
                      'setup.py.template',
                      ],
        },
        zip_safe=False,
        include_package_data=True,
        exclude_package_data={'': ['README.txt'],
                              # 'amici': ['src/*']
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
