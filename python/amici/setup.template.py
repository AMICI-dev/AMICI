"""AMICI model package setup"""

from setuptools import find_packages
from distutils.core import setup, Extension
from distutils import sysconfig
import os
from amici import amici_path, hdf5_enabled

from amici.setuptools import (getBlasConfig,
                              getHdf5Config,
                              addCoverageFlagsIfRequired,
                              addDebugFlagsIfRequired)

os.chdir(os.path.dirname(os.path.abspath(__file__)))

def getModelSources():
    """Get list of source files for the amici base library"""
    import glob
    modelSources = glob.glob('*.cpp')
    try:
        modelSources.remove('main.cpp')
    except ValueError:
        pass
    return modelSources


def getAmiciLibs():
    """Get list of libraries for the amici base library"""
    return ['amici',
            'sundials', 'suitesparse'
            #'sundials_nvecserial', 'sundials_cvodes', 'sundials_idas',
            #'klu', 'colamd', 'btf', 'amd', 'suitesparseconfig'
            ]

cxx_flags = ['-std=c++11']
linker_flags = []

addCoverageFlagsIfRequired(cxx_flags, linker_flags)
addDebugFlagsIfRequired(cxx_flags, linker_flags)

h5pkgcfg = getHdf5Config()

blaspkgcfg = getBlasConfig()
linker_flags.extend(blaspkgcfg.get('extra_link_args', []))

libraries = [*getAmiciLibs(), *blaspkgcfg['libraries']]
if hdf5_enabled:
    libraries.extend(['hdf5_hl_cpp', 'hdf5_hl', 'hdf5_cpp', 'hdf5'])

sources = ['swig/TPL_MODELNAME.i', *getModelSources()]


# Remove the "-Wstrict-prototypes" compiler option, which isn't valid for
# C++ to fix warnings.
cfg_vars = sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if isinstance(value, str):
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")

# compiler and linker flags for libamici
if 'AMICI_CXXFLAGS' in os.environ:
    cxx_flags.extend(os.environ['AMICI_CXXFLAGS'].split(' '))
if 'AMICI_LDFLAGS' in os.environ:
    linker_flags.extend(os.environ['AMICI_LDFLAGS'].split(' '))

# Build shared object
model_module = Extension('TPL_MODELNAME._TPL_MODELNAME',
                         sources=sources,
                         include_dirs=[os.getcwd(),
                                       os.path.join(amici_path, 'include'),
                                       os.path.join(amici_path, 'ThirdParty/sundials/include'),
                                       os.path.join(amici_path, 'ThirdParty/SuiteSparse/include'),
                                       *h5pkgcfg['include_dirs'],
                                       *blaspkgcfg['include_dirs']
                                       ],
                         libraries = libraries,
                         library_dirs=[
                             *h5pkgcfg['library_dirs'],
                             *blaspkgcfg['library_dirs'],
                             os.path.join(amici_path, 'libs')],
                         swig_opts=['-c++', '-modern', '-outdir', 'TPL_MODELNAME',
                                    '-I%s' % os.path.join(amici_path, 'swig'),
                                    '-I%s' % os.path.join(amici_path, 'include'),],
                         extra_compile_args=cxx_flags,
                         extra_link_args=linker_flags
                         )

# Install
setup(
    name='TPL_MODELNAME',
    version='TPL_PACKAGE_VERSION',
    description='AMICI-generated module for model TPL_MODELNAME',
    url='https://github.com/ICB-DCM/AMICI',
    author='model-author-todo',
    author_email='model-author-todo',
    #license = 'BSD',
    ext_modules=[model_module],
    packages=find_packages(),
    install_requires=['amici==TPL_AMICI_VERSION'],
    extras_require={'wurlitzer': ['wurlitzer']},
    python_requires='>=3.6',
    package_data={
    },
    zip_safe = False,
    include_package_data=True,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        #'License :: OSI Approved :: BSD License',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS :: MacOS X',
        'Programming Language :: Python',
        'Programming Language :: C++',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
