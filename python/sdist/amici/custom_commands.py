"""Custom setuptools commands for AMICI installation"""

import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import cmake_build_extension
from setuptools.command.build_ext import build_ext
from setuptools.command.build_py import build_py
from setuptools.command.develop import develop
from setuptools.command.install import install
from setuptools.command.install_lib import install_lib
from setuptools.command.sdist import sdist

from amici.swig import fix_typehints

# typehints
Library = Tuple[str, Dict[str, List[str]]]


class AmiciInstall(install):
    """Custom install to handle extra arguments"""

    print("running AmiciInstall")

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


class AmiciInstallLib(install_lib):
    """Custom install to allow preserving of debug symbols"""

    def run(self):
        """strip debug symbols

        Returns:

        """
        print("running AmiciInstallLib")

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


# TODO REMOVEME
class AmiciBuildExt(build_ext):
    """Custom build_ext to allow keeping otherwise temporary static libs"""

    def run(self):
        """Copy the generated clibs to the extensions folder to be included in
        the wheel
        """

        print("running AmiciBuildExt")

        no_clibs = 'develop' in self.distribution.command_obj \
                   and self.get_finalized_command('develop').no_clibs
        no_clibs |= 'install' in self.distribution.command_obj \
                    and self.get_finalized_command('install').no_clibs

        if no_clibs:
            # Nothing to build
            return

        if not self.dry_run and self.distribution.has_c_libraries():
            # Module build directory where we want to copy the generated
            # libs to
            if self.inplace == 0:
                build_dir = self.build_lib
            else:
                build_dir = os.getcwd()
            target_dir = os.path.join(build_dir, 'amici', 'libs')
            self.mkpath(target_dir)

            swig_outdir = os.path.join(os.path.abspath(build_dir), "amici")
            swig_py_module_path = os.path.join(swig_outdir, 'amici.py')
            print("updating typehints")
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

        print("running AmiciSDist")

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
            print(e)


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
                print(f"Changed {field} for {lib[0]} with {compiler_type} "
                         f"to {lib[1][field]}")
            except KeyError:
                # No compiler-specific options set
                pass


class AmiciBuildPy(build_py):
    def run(self):
        # We need build_ext before build_py, that all artifacts will be
        # copied from the build dir
        self.run_command("build_clib")
        self.run_command("build_ext")
        return super().run()



class AmiciBuildCMakeExtension(cmake_build_extension.BuildExtension):
    def build_extension(self, ext: cmake_build_extension.CMakeExtension) -> None:
        print("-" * 20, ext.name, "-" * 20, file=sys.stderr)

        if self.inplace == 0:
            build_dir = self.build_lib
        else:
            build_dir = os.getcwd()
        print("*" * 50)
        print(sys.implementation)
        print("*" * 50)
        for x in dir(self):
            print(x, getattr(self, x))
        print("*" * 50)

        import sysconfig
        print(sysconfig.get_config_vars())
        if sys.platform == "win32":
            # account for e.g. build/lib.win-amd64-cpython-38
            build_dir = Path(build_dir, "build", f"lib.{sysconfig.get_platform()}-{sys.implementation.cache_tag}").absolute().as_posix()
        else:
            build_dir = str(Path(build_dir).absolute())

        #ext.cmake_configure_options = [x.replace("${suitesparse_root}", clib_dir) for x in ext.cmake_configure_options]
        ext.cmake_configure_options = [
            x.replace("${build_dir}", build_dir) for x in
            ext.cmake_configure_options]
        cmake_build_extension.BuildExtension.build_extension(self, ext)

        print("-" * 40, file=sys.stderr)
