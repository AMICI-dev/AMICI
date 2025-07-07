"""Functions related to SWIG or SWIG-generated code"""

from __future__ import annotations
import ast
import contextlib
import re


class TypeHintFixer(ast.NodeTransformer):
    """Replaces SWIG-generated C++ typehints by corresponding Python types."""

    mapping = {
        "void": None,
        "double": ast.Name("float"),
        "int": ast.Name("int"),
        "long": ast.Name("int"),
        "ptrdiff_t": ast.Name("int"),
        "size_t": ast.Name("int"),
        "bool": ast.Name("bool"),
        "boolean": ast.Name("bool"),
        "std::unique_ptr< amici::Solver >": ast.Constant("Solver"),
        "amici::InternalSensitivityMethod": ast.Constant(
            "InternalSensitivityMethod"
        ),
        "amici::InterpolationType": ast.Constant("InterpolationType"),
        "amici::LinearMultistepMethod": ast.Constant("LinearMultistepMethod"),
        "amici::LinearSolver": ast.Constant("LinearSolver"),
        "amici::Model *": ast.Constant("Model"),
        "amici::Model const *": ast.Constant("Model"),
        "amici::NewtonDampingFactorMode": ast.Constant(
            "NewtonDampingFactorMode"
        ),
        "amici::NonlinearSolverIteration": ast.Constant(
            "NonlinearSolverIteration"
        ),
        "amici::ObservableScaling": ast.Constant("ObservableScaling"),
        "amici::ParameterScaling": ast.Constant("ParameterScaling"),
        "amici::RDataReporting": ast.Constant("RDataReporting"),
        "amici::SensitivityMethod": ast.Constant("SensitivityMethod"),
        "amici::SensitivityOrder": ast.Constant("SensitivityOrder"),
        "amici::Solver *": ast.Constant("Solver"),
        "amici::SteadyStateSensitivityMode": ast.Constant(
            "SteadyStateSensitivityMode"
        ),
        "amici::realtype": ast.Name("float"),
        "DoubleVector": ast.Name("Sequence[float]"),
        "BoolVector": ast.Name("Sequence[bool]"),
        "IntVector": ast.Name("Sequence[int]"),
        "StringVector": ast.Name("Sequence[str]"),
        "std::string": ast.Name("str"),
        "std::string const &": ast.Name("str"),
        "std::unique_ptr< amici::ExpData >": ast.Constant("ExpData"),
        "std::unique_ptr< amici::ReturnData >": ast.Constant("ReturnData"),
        "std::vector< amici::ParameterScaling,"
        "std::allocator< amici::ParameterScaling > > const &": ast.Constant(
            "ParameterScalingVector"
        ),
        "H5::H5File": None,
    }

    def visit_FunctionDef(self, node):
        # convert type/rtype from docstring to annotation, if possible.
        #  those may be c++ types, not valid in python, that need to be
        #  converted to python types below.
        self._annotation_from_docstring(node)

        # Has a return type annotation?
        if node.returns and isinstance(node.returns, ast.Constant):
            node.returns = self._new_annot(node.returns.value)

        # Has arguments?
        if node.args.args:
            for arg in node.args.args:
                if not arg.annotation:
                    continue
                if not isinstance(arg.annotation, ast.Constant):
                    # there is already proper annotation
                    continue

                arg.annotation = self._new_annot(arg.annotation.value)
        return node

    def _new_annot(self, old_annot: str | ast.Name):
        if isinstance(old_annot, ast.Name):
            old_annot = old_annot.id

        with contextlib.suppress(KeyError):
            return self.mapping[old_annot]

        # std::vector size type
        if re.match(r"std::vector< .* >::(?:size|difference)_type", old_annot):
            return ast.Name("int")

        # std::vector value type
        if (
            value_type := re.sub(
                r"std::vector< (.*) >::value_type(?: const &)?",
                r"\1",
                old_annot,
            )
        ) in self.mapping:
            return self.mapping[value_type]

        # std::vector
        if (
            value_type := re.sub(
                r"std::vector< (.*),std::allocator< \1 > >(?: const &)?",
                r"\1",
                old_annot,
            )
        ) in self.mapping:
            value_type_annot = self.mapping[value_type]
            if isinstance(value_type_annot, ast.Constant):
                return ast.Name(f"tuple['{value_type_annot.value}']")
            if isinstance(value_type_annot, ast.Name):
                return ast.Name(f"tuple[{value_type_annot.id}]")

        return ast.Constant(old_annot)

    def _annotation_from_docstring(self, node: ast.FunctionDef):
        """Add annotations based on docstring.

        If any argument or return type of the function is not annotated, but
        the corresponding docstring contains a type hint (``:rtype:`` or
        ``:type:``), the type hint is used as the annotation.

        Swig sometimes generates ``:type solver: :py:class:`Solver`` instead of
        ``:type solver: Solver``. Those need special treatment.

        Overloaded functions are skipped.
        """
        docstring = ast.get_docstring(node, clean=False)
        if not docstring or "*Overload 1:*" in docstring:
            # skip overloaded methods
            return

        docstring = docstring.split("\n")
        lines_to_remove = set()

        for line_no, line in enumerate(docstring):
            if type_str := self.extract_rtype(line):
                # handle `:rtype:`
                node.returns = ast.Constant(type_str)
                lines_to_remove.add(line_no)
                continue

            arg_name, type_str = self.extract_type(line)
            if arg_name is not None:
                # handle `:type ...:`
                for arg in node.args.args:
                    if arg.arg == arg_name:
                        arg.annotation = ast.Constant(type_str)
                        lines_to_remove.add(line_no)

        if lines_to_remove:
            # Update docstring with type annotations removed
            assert isinstance(node.body[0].value, ast.Constant)
            new_docstring = "\n".join(
                line
                for line_no, line in enumerate(docstring)
                if line_no not in lines_to_remove
            )
            node.body[0].value = ast.Str(new_docstring)

    @staticmethod
    def extract_type(line: str) -> tuple[str, str] | tuple[None, None]:
        """Extract argument name and type string from ``:type:`` docstring
        line."""
        match = re.match(r"\s*:type\s+(\w+):\s+(.+?)(?:, optional)?\s*$", line)
        if not match:
            return None, None

        arg_name = match.group(1)

        # get rid of any :py:class`...` in the type string if necessary
        if not match.group(2).startswith(":py:"):
            return arg_name, match.group(2)

        match = re.match(r":py:\w+:`(.+)`", match.group(2))
        assert match
        return arg_name, match.group(1)

    @staticmethod
    def extract_rtype(line: str) -> str | None:
        """Extract type string from ``:rtype:`` docstring line."""
        match = re.match(r"\s*:rtype:\s+(.+)\s*$", line)
        if not match:
            return None

        # get rid of any :py:class`...` in the type string if necessary
        if not match.group(1).startswith(":py:"):
            return match.group(1)

        match = re.match(r":py:\w+:`(.+)`", match.group(1))
        assert match
        return match.group(1)


def fix_typehints(infilename, outfilename):
    """Change SWIG-generated C++ typehints to Python typehints"""
    # file -> AST
    with open(infilename) as f:
        source = f.read()
    parsed_source = ast.parse(source)

    # Change AST
    fixer = TypeHintFixer()
    parsed_source = fixer.visit(parsed_source)

    # AST -> file
    with open(outfilename, "w") as f:
        f.write(ast.unparse(parsed_source))
