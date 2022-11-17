"""Functions for downloading/building/finding SWIG"""
import ast
import contextlib
import os
import re
import subprocess
from typing import Tuple


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
        'double': ast.Name('float'),
        'int': ast.Name('int'),
        'long': ast.Name('int'),
        'ptrdiff_t': ast.Name('int'),
        'size_t': ast.Name('int'),
        'bool': ast.Name('bool'),
        'std::unique_ptr< amici::Solver >': ast.Constant('Solver'),
        'amici::InternalSensitivityMethod':
            ast.Constant('InternalSensitivityMethod'),
        'amici::InterpolationType': ast.Constant('InterpolationType'),
        'amici::LinearMultistepMethod': ast.Constant('LinearMultistepMethod'),
        'amici::LinearSolver': ast.Constant('LinearSolver'),
        'amici::Model *': ast.Constant('Model'),
        'amici::Model const *': ast.Constant('Model'),
        'amici::NewtonDampingFactorMode':
            ast.Constant('NewtonDampingFactorMode'),
        'amici::NonlinearSolverIteration':
            ast.Constant('NonlinearSolverIteration'),
        'amici::ObservableScaling': ast.Constant('ObservableScaling'),
        'amici::ParameterScaling': ast.Constant('ParameterScaling'),
        'amici::RDataReporting': ast.Constant('RDataReporting'),
        'amici::SensitivityMethod': ast.Constant('SensitivityMethod'),
        'amici::SensitivityOrder': ast.Constant('SensitivityOrder'),
        'amici::Solver *': ast.Constant('Solver'),
        'amici::SteadyStateSensitivityMode':
            ast.Constant('SteadyStateSensitivityMode'),
        'amici::realtype': ast.Name('float'),
        'DoubleVector': ast.Constant('Sequence[float]'),
        'IntVector': ast.Name('Sequence[int]'),
        'std::string': ast.Name('str'),
        'std::string const &': ast.Name('str'),
        'std::unique_ptr< amici::ExpData >': ast.Constant('ExpData'),
        'std::unique_ptr< amici::ReturnData >': ast.Constant('ReturnData'),
        'std::vector< amici::ParameterScaling,'
        'std::allocator< amici::ParameterScaling > > const &':
            ast.Constant('ParameterScalingVector')
    }

    def visit_FunctionDef(self, node):
        # Has a return type annotation?
        if node.returns:
            node.returns = self._new_annot(node.returns.value)

        # Has arguments?
        if node.args.args:
            for arg in node.args.args:
                if not arg.annotation:
                    continue
                arg.annotation = self._new_annot(arg.annotation.value)
        return node

    def _new_annot(self, old_annot: str):
        with contextlib.suppress(KeyError):
            return self.mapping[old_annot]

        # std::vector size type
        if re.match(r"std::vector< .* >::(?:size|difference)_type", old_annot):
            return ast.Name("int")

        # std::vector value type
        if (value_type := re.sub(
                r'std::vector< (.*) >::value_type(?: const &)?',
                r'\1', old_annot)) in self.mapping:
            return self.mapping[value_type]

        # std::vector
        if (value_type := re.sub(
                r'std::vector< (.*),std::allocator< \1 > >(?: const &)?',
                r'\1', old_annot)) in self.mapping:
            value_type_annot = self.mapping[value_type]
            if isinstance(value_type_annot, ast.Constant):
                return ast.Name(f"Tuple['{value_type_annot.value}']")
            if isinstance(value_type_annot, ast.Name):
                return ast.Name(f"Tuple[{value_type_annot.id}]")

        return ast.Constant(old_annot)


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
