"""Functions for downloading/building/finding SWIG"""

from typing import Tuple
import ast
import os
import subprocess
import re


def find_swig() -> str:
    """Get name and version of SWIG executable

    We need version >=3.0. Probably we should try some default paths and names,
    but this should do the trick for now.

    Debian/Ubuntu systems have swig3.0 ('swig' is older versions),
    OSX has swig 3.0 as 'swig'.
    """

    candidates = ['swig4.0', 'swig3.0', 'swig']
    # Environment variable has priority
    if 'SWIG' in os.environ:
        candidates.insert(0, os.environ['SWIG'])

    for candidate in candidates:
        if swig_works(candidate):
            return candidate

    raise RuntimeError(
        "Unable to find SWIG executable with default names. "
        "Ensure you have SWIG installed, e.g. by "
        "`sudo apt install swig` or `brew install swig`. "
        "As non-root user, you can install SWIG using "
        "https://github.com/AMICI-dev/AMICI/blob/master/scripts/"
        "downloadAndBuildSwig.sh, or by following the "
        "instructions at http://www.swig.org/Doc4.0/"
        "SWIGDocumentation.html#Preface_installation. "
        "If was not found despite being installed, set the SWIG"
        " environment variable to the full path of the correct "
        "executable."
    )


def swig_works(swig: str, verbose: bool = True) -> bool:
    """Test if `swig` looks like a working SWIG executable."""

    try:
        # For python3.6 compatibility no `capture_output=True`
        result = subprocess.run([swig, '-version'],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    except (FileNotFoundError, PermissionError):
        if verbose:
            print(f'Testing SWIG executable {swig}... FAILED.')
        return False

    if verbose:
        if result.returncode == 0:
            print(f'Testing SWIG executable {swig}... SUCCEEDED.')
        else:
            print(f'Testing SWIG executable {swig}... FAILED.')

    return result.returncode == 0


def get_swig_version(swig_exe: str) -> Tuple:
    """Determine version of the given SWIG executable

    Returns:
        Version tuple
    """
    result = subprocess.run([swig_exe, '-version'],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    assert result.returncode == 0

    version = re.sub(r'(?s).*Version\s+([\S]+).*', r'\1',
                     result.stdout.decode('utf-8'))

    return tuple(int(x) for x in version.split('.'))


class TypeHintFixer(ast.NodeTransformer):
    """Replaces SWIG-generated C++ typehints by corresponding Python types"""

    mapping = {
        'void': None,
        'std::unique_ptr< amici::Solver >': 'amici.Solver',
        'amici::InternalSensitivityMethod': 'amici.InternalSensitivityMethod',
        'amici::InterpolationType': 'amici.InterpolationType',
        'amici::LinearMultistepMethod': 'amici.LinearMultistepMethod',
        'amici::LinearSolver': 'amici.LinearSolver',
        'amici::Model *': 'amici.Model',
        'amici::Model const *': 'amici.Model',
        'amici::NewtonDampingFactorMode': 'amici.NewtonDampingFactorMode',
        'amici::NonlinearSolverIteration': 'amici.NonlinearSolverIteration',
        'amici::RDataReporting': 'amici.RDataReporting',
        'amici::SensitivityMethod': 'amici.SensitivityMethod',
        'amici::SensitivityOrder': 'amici.SensitivityOrder',
        'amici::Solver *': 'amici.Solver',
        'amici::SteadyStateSensitivityMode': 'amici.SteadyStateSensitivityMode',
        'amici::realtype': 'float',
        'DoubleVector': 'numpy.ndarray',
        'IntVector': 'List[int]',
        'std::string': 'str',
        'std::string const &': 'str',
        'std::unique_ptr< amici::ExpData >': 'amici.ExpData',
        'std::unique_ptr< amici::ReturnData >': 'amici.ReturnData',
    }

    def visit_FunctionDef(self, node):
        # Has a return type annotation?
        if node.returns:
            node.returns.value = self._new_annot(node.returns.value)

        # Has arguments?
        if node.args.args:
            for arg in node.args.args:
                if not arg.annotation:
                    continue
                arg.annotation.value = self._new_annot(arg.annotation.value)
        return node

    def _new_annot(self, old_annot):
        return self.mapping.get(old_annot, old_annot)


def fix_typehints(infilename, outfilename):
    """Change SWIG-generated C++ typehints to Python typehints"""
    # Only available from Python3.9
    if not getattr(ast, 'unparse', None):
        return

    # file -> AST
    with open(infilename, 'r') as f:
        source = f.read()
    parsed_source = ast.parse(source)

    # Change AST
    fixer = TypeHintFixer()
    parsed_source = fixer.visit(parsed_source)

    # AST -> file
    with open(outfilename, 'w') as f:
        f.write(ast.unparse(parsed_source))
