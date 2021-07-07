"""Custom setuptools commands for AMICI installation"""

import glob
import os
import subprocess
import sys
from shutil import copyfile
from typing import Dict, List, Tuple

from amici.swig import fix_typehints
from amici.setuptools import generate_swig_interface_files
from setuptools.command.build_clib import build_clib
from setuptools.command.build_ext import build_ext
from setuptools.command.develop import develop
from setuptools.command.install import install
from setuptools.command.install_lib import install_lib
from setuptools.command.sdist import sdist
from distutils import log

# typehints
Library = Tuple[str, Dict[str, List[str]]]


class AmiciInstall(install):
    """Custom install to handle extra arguments"""

    log.debug("running AmiciInstall")

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


class AmiciBuildCLib(build_clib):
    """Custom build_clib"""

    def run(self):
        log.debug("running AmiciBuildCLib")

        # Always force recompilation. The way setuptools/distutils check for
        # whether sources require recompilation is not reliable and may lead
        # to crashes or wrong results. We rather compile once too often...
        self.force = True

        build_clib.run(self)

    def build_libraries(self, libraries: List[Library]):
        log.debug("running AmiciBuildCLib.build_libraries")

        no_clibs = 'develop' in self.distribution.command_obj \
                   and self.get_finalized_command('develop').no_clibs
        no_clibs |= 'install' in self.distribution.command_obj \
                    and self.get_finalized_command('install').no_clibs

        if no_clibs:
            return

        # Override for parallel compilation
        import distutils.ccompiler
        distutils.ccompiler.CCompiler.compile = compile_parallel

        # Work-around for compiler-specific build options
        set_compiler_specific_library_options(
            libraries, self.compiler.compiler_type)

        # Monkey-patch setuptools, to force recompilation of library sources
        # --force does not work as expected

        # need full import here, not module-level imported build_clib
        import setuptools.command.build_clib
        # the patched function may return anything but `([], [])` to trigger
        # recompilation
        setuptools.command.build_clib.newer_pairwise_group = lambda *_: None

        build_clib.build_libraries(self, libraries)


class AmiciDevelop(develop):
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
        log.debug("running AmiciDevelop")

        if not self.no_clibs:
            self.get_finalized_command('build_clib').run()

        develop.run(self)


class AmiciInstallLib(install_lib):
    """Custom install to allow preserving of debug symbols"""

    def run(self):
        """strip debug symbols

        Returns:

        """
        log.debug("running AmiciInstallLib")

        if 'ENABLE_AMICI_DEBUGGING' in os.environ \
                and os.environ['ENABLE_AMICI_DEBUGGING'] == 'TRUE' \
                and sys.platform == 'darwin':
            search_dir = os.path.join(os.getcwd(), self.build_dir, 'amici')
            for file in os.listdir(search_dir):
                if file.endswith('.so'):
                    subprocess.run(['dsymutil', os.path.join(search_dir, file),
                                    '-o',
                                    os.path.join(search_dir, file + '.dSYM')])

        # Continue with the actual installation
        install_lib.run(self)


class AmiciBuildExt(build_ext):
    """Custom build_ext to allow keeping otherwise temporary static libs"""

    def build_extension(self, ext):
        # Work-around for compiler-specific build options
        set_compiler_specific_extension_options(
            ext, self.compiler.compiler_type)

        build_ext.build_extension(self, ext)

    def run(self):
        """Copy the generated clibs to the extensions folder to be included in
        the wheel
        """

        log.debug("running AmiciBuildExt")

        no_clibs = 'develop' in self.distribution.command_obj \
                   and self.get_finalized_command('develop').no_clibs
        no_clibs |= 'install' in self.distribution.command_obj \
                    and self.get_finalized_command('install').no_clibs

        if no_clibs:
            # Nothing to build
            return

        if not self.dry_run and self.distribution.has_c_libraries():
            # get the previously built static libraries
            build_clib = self.get_finalized_command('build_clib')
            libraries = build_clib.get_library_names() or []

            # Module build directory where we want to copy the generated
            # libs to
            if self.inplace == 0:
                build_dir = self.build_lib
            else:
                build_dir = os.getcwd()
            target_dir = os.path.join(build_dir, 'amici', 'libs')
            self.mkpath(target_dir)

            # Copy the generated libs
            for lib in libraries:
                libfilenames = glob.glob(
                    f"{build_clib.build_clib}{os.sep}*{lib}.*")
                assert len(libfilenames) == 1, \
                    f"Found unexpected number of files: {libfilenames}"
                src = libfilenames[0]
                dest = os.path.join(target_dir, os.path.basename(src))
                log.info(f"copying {src} -> {dest}")
                copyfile(src, dest)

            swig_outdir = os.path.join(os.path.abspath(build_dir), "amici")
            generate_swig_interface_files(swig_outdir=swig_outdir)
            swig_py_module_path = os.path.join(swig_outdir, 'amici.py')
            log.debug("updating typehints")
            fix_typehints(swig_py_module_path, swig_py_module_path)

        # Always force recompilation. The way setuptools/distutils check for
        # whether sources require recompilation is not reliable and may lead
        # to crashes or wrong results. We rather compile once too often...
        self.force = True

        # Continue with the actual extension building
        build_ext.run(self)


class AmiciSDist(sdist):
    """Customized creation of source distribution"""

    def run(self):
        """Setuptools entry-point"""

        log.debug("running AmiciSDist")

        save_git_version()

        sdist.run(self)


def save_git_version():
    """Create file with extended version string

    This requires git. We assume that whoever creates the sdist will work
    inside a valid git repository.

    Returns:

    """
    with open(os.path.join("amici", "git_version.txt"), "w") as f:
        try:
            cmd = ['git', 'describe', '--abbrev=4', '--dirty=-dirty',
                   '--always', '--tags']
            subprocess.run(cmd, stdout=f)
        except Exception as e:
            log.warn(e)


def set_compiler_specific_library_options(
        libraries: List[Library],
        compiler_type: str) -> None:
    """Set compiler-specific library options.

    C/C++-libraries for setuptools/distutils are provided as dict containing
    entries for 'sources', 'macros', 'cflags', etc.
    As we don't know the compiler type at the stage of calling
    ``setuptools.setup`` and as there is no other apparent way to set
    compiler-specific options, we elsewhere extend the dict with additional
    fields ${original_field}_${compiler_class}, and add the additional
    compiler-specific options here, at a stage when the compiler has been
    determined by distutils.

    Arguments:
        libraries:
            List of libraries as passed as ``libraries`` argument to
            ``setuptools.setup`` and ``setuptools.build_ext.build_extension``.
            This is modified in place.
        compiler_type:
            Compiler type, as defined in
            ``distutils.ccompiler.compiler.compiler_class``, (e.g. 'unix',
            'msvc', 'mingw32').
    """

    for lib in libraries:
        for field in ['cflags', 'sources', 'macros']:
            try:
                lib[1][field] += lib[1][f'{field}_{compiler_type}']
                log.info(f"Changed {field} for {lib[0]} with {compiler_type} "
                         f"to {lib[1][field]}")
            except KeyError:
                # No compiler-specific options set
                pass


def set_compiler_specific_extension_options(
        ext: 'setuptools.Extension',
        compiler_type: str) -> None:
    """Set compiler-specific extension build options.

    Same game as in ``set_compiler_specific_library_options``, except that
    here we look for compiler-specific class attributes.

    Arguments:
        ext: setuptools/distutils extension object
        compiler_type: Compiler type
    """
    for attr in ['extra_compile_args', 'extra_link_args']:
        try:
            new_value = getattr(ext, attr) + \
                        getattr(ext, f'{attr}_{compiler_type}')
            setattr(ext, attr, new_value)
            log.info(f"Changed {attr} for {compiler_type} to {new_value}")
        except AttributeError:
            # No compiler-specific options set
            pass

