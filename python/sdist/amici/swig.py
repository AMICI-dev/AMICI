"""Functions related to SWIG or SWIG-generated code"""
import ast
import contextlib
import re


class TypeHintFixer(ast.NodeTransformer):
    """Replaces SWIG-generated C++ typehints by corresponding Python types"""

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
    }

    def visit_FunctionDef(self, node):
        self._annotation_from_docstring(node)

        # Has a return type annotation?
        if node.returns and isinstance(node.returns, ast.Constant):
            node.returns = self._new_annot(node.returns.value)

        # Has arguments?
        if node.args.args:
            for arg in node.args.args:
                if not arg.annotation:
                    continue
                if isinstance(arg.annotation, ast.Name):
                    # there is already proper annotation
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
        """
        docstring = ast.get_docstring(node, clean=False)
        if not docstring or "*Overload 1:*" in docstring:
            # skip overloaded methods
            return

        docstring = docstring.split("\n")
        lines_to_remove = set()

        for line_no, line in enumerate(docstring):
            if (
                match := re.match(
                    r"\s*:rtype:\s*(?::py:class:`)?(\w+)`?\s+$", line
                )
            ) and not match.group(1).startswith(":"):
                node.returns = ast.Constant(match.group(1))
                lines_to_remove.add(line_no)

            if (
                match := re.match(
                    r"\s*:type\s*(\w+):\W*(?::py:class:`)?(\w+)`?\s+$", line
                )
            ) and not match.group(1).startswith(":"):
                for arg in node.args.args:
                    if arg.arg == match.group(1):
                        arg.annotation = ast.Constant(match.group(2))
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


def fix_typehints(infilename, outfilename):
    """Change SWIG-generated C++ typehints to Python typehints"""
    # Only available from Python3.9
    if not getattr(ast, "unparse", None):
        return

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
