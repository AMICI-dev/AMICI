"""Custom setuptools commands for AMICI installation"""

import glob
import os
import sys
import subprocess
from shutil import copyfile

from setuptools.command.build_ext import build_ext
from setuptools.command.sdist import sdist
from setuptools.command.install_lib import install_lib
from setuptools.command.develop import develop
from setuptools.command.install import install
from setuptools.command.build_clib import build_clib

from amici.setuptools import generateSwigInterfaceFiles


class my_install(install):
    """Custom install to handle extra arguments"""

    # Passing --no-clibs allows to install the Python-only part of AMICI
    user_options = install.user_options + [
        ('no-clibs', None, "Don't build AMICI C++ extension"),
    ]

    def initialize_options(self):
        install.initialize_options(self)
        self.no_clibs = False

    def finalize_options(self):
        if self.no_clibs:
            self.no_clibs = True
        install.finalize_options(self)

    def run(self):
        generateSwigInterfaceFiles()
        install.run(self)


def compile_parallel(self, sources, output_dir=None, macros=None,
            include_dirs=None, debug=0, extra_preargs=None,
            extra_postargs=None, depends=None):
    """Parallelized version of distutils.ccompiler.compile"""
    macros, objects, extra_postargs, pp_opts, build = \
        self._setup_compile(output_dir, macros, include_dirs, sources,
                            depends, extra_postargs)
    cc_args = self._get_cc_args(pp_opts, debug, extra_preargs)

    # parallel compilation
    num_threads = 1
    if 'AMICI_PARALLEL_COMPILE' in os.environ:
        max_threads = int(os.environ['AMICI_PARALLEL_COMPILE'])
        num_threads = min(len(objects), max_threads)
        num_threads = max(1, num_threads)

    def _single_compile(obj):
        try:
            src, ext = build[obj]
        except KeyError:
            return
        self._compile(obj, src, ext, cc_args, extra_postargs, pp_opts)

    if num_threads > 1:
        import multiprocessing.pool
        # convert to list, imap is evaluated on-demand
        list(multiprocessing.pool.ThreadPool(num_threads).imap(
            _single_compile, objects))
    else:
        for obj in objects:
            _single_compile(obj)

    return objects


class my_build_clib(build_clib):
    """Custom build_clib"""

    def run(self):
        # Always force recompilation. The way setuptools/distutils check for
        # whether sources require recompilation is not reliable and may lead
        # to crashes or wrong results. We rather compile once too often...
        self.force = True

        build_clib.run(self)


    def build_libraries(self, libraries):
        no_clibs = self.get_finalized_command('develop').no_clibs
        no_clibs |= self.get_finalized_command('install').no_clibs

        if no_clibs:
            return

        # Override for parallel compilation
        import distutils.ccompiler
        distutils.ccompiler.CCompiler.compile = compile_parallel

        # start new code
        compilerType = self.compiler.compiler_type
        for lib in libraries:
            try:
               lib[1]['cflags'] = lib[1]['cflags'] + lib[1]['cflags_'+compilerType]
            except KeyError: None
        # end new code

        build_clib.build_libraries(self, libraries)

class my_develop(develop):
    """Custom develop to build clibs"""

    # Passing --no-clibs allows to install the Python-only part of AMICI
    user_options = develop.user_options + [
        ('no-clibs', None, "Don't build AMICI C++ extension"),
    ]

    def initialize_options(self):
        develop.initialize_options(self)
        self.no_clibs = False

    def finalize_options(self):
        if self.no_clibs:
            self.no_clibs = True
        develop.finalize_options(self)

    def run(self):
        if not self.no_clibs:
            generateSwigInterfaceFiles()
            self.run_command('build')

        develop.run(self)


class my_install_lib(install_lib):
    """Custom install to allow preserving of debug symbols"""
    def run(self):
        """strip debug symbols

        Returns:

        """
        if 'ENABLE_AMICI_DEBUGGING' in os.environ \
                and os.environ['ENABLE_AMICI_DEBUGGING'] == 'TRUE' \
                and sys.platform == 'darwin':
            search_dir = os.path.join(os.getcwd(),self.build_dir,'amici')
            for file in os.listdir(search_dir):
                if file.endswith('.so'):
                    subprocess.run(['dsymutil',os.path.join(search_dir,file),
                                    '-o',os.path.join(search_dir,file + '.dSYM')])


        # Continue with the actual installation
        install_lib.run(self)


class my_build_ext(build_ext):
    """Custom build_ext to allow keeping otherwise temporary static libs"""

    def build_extension(self, ext):
        # do replacements here
        compilerType = self.compiler.compiler_type
        # this requires more work, or maybe deletation of the commented code
        #try: ext.extra_compile_args = ext.extra_compile_args, ext.(extra_compile_args_+compilerType)]
        #except KeyError: None
        if compilerType == 'msvc':
            ECA = ext.extra_compile_args
            for index in range(len(ECA)):
                ECA[index] = ECA[index].replace('-std=c++14', '/std:c++14')
            ext.extra_compile_args = ECA

        build_ext.build_extension(self, ext)

    def run(self):
        """Copy the generated clibs to the extensions folder to be included in the wheel

        Returns:

        """
        no_clibs = self.get_finalized_command('develop').no_clibs
        no_clibs |= self.get_finalized_command('install').no_clibs

        if no_clibs:
            return

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
                assert len(libfilenames) == 1, \
                    "Found unexpected number of files: " % libfilenames

                copyfile(libfilenames[0],
                         os.path.join(target_dir, os.path.basename(libfilenames[0])))

        # Always force recompilation. The way setuptools/distutils check for
        # whether sources require recompilation is not reliable and may lead
        # to crashes or wrong results. We rather compile once too often...
        self.force = True

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
