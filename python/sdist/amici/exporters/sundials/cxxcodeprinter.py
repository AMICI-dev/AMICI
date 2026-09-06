"""C++ code generation"""

import itertools
import os
import re
import warnings
from collections.abc import Iterable, Sequence

import sympy as sp
from sympy.codegen.rewriting import Optimization, optimize
from sympy.printing.cxx import CXX11CodePrinter
from sympy.utilities.iterables import numbered_symbols
from toposort import toposort

from amici.importers.utils import amici_time_symbol, symbol_with_assumptions


def _mangle(name: str) -> str:
    """Make a model-derived identifier safe as a C++ local variable name.

    Appends `_` so it can never equal a real keyword/macro. Collapses any
    SBML-legal `__` run first, and uses `v` instead of `_` as the marker for
    names already ending in `_`, so the result never contains `__` either.
    This is not injective and does not guarantee a unique result on its own
    (e.g. `"a__b"` and `"a_b"` both collapse to the same string) -- callers
    that need uniqueness across a whole model handle that separately.
    """
    name = re.sub(r"_{2,}", "_", name)
    return f"{name}v" if name.endswith("_") else f"{name}_"


class AmiciCxxCodePrinter(CXX11CodePrinter):
    """
    C++ code printer

    Attributes
    ----------
    extract_cse:
        Whether to extract common subexpression during code printing.
        Currently controlled by environment variable ``AMICI_EXTRACT_CSE``.
    optimizations:
        Iterable of :class:`sympy.codegen.rewriting.Optimization`s to optimize
        generated code (e.g. :data:`sympy.codegen.rewriting.Optimization` for
        optimizations, such as ``log(1 + x)`` --> ``logp1(x)``).
        Applying these optimizations is potentially quite costly.
    """

    optimizations: Iterable[Optimization] = ()

    def __init__(self, extract_cse: bool | None = None):
        """Create code printer"""
        super().__init__()

        # extract common subexpressions in matrix functions?
        self.extract_cse = (
            os.getenv("AMICI_EXTRACT_CSE", "0").lower()
            in (
                "1",
                "on",
                "true",
            )
            if extract_cse is None
            else extract_cse
        )

        # Floating-point optimizations
        # e.g., log(1 + x) --> logp1(x)
        if self.optimizations:
            self._fpoptimizer = lambda x: optimize(x, self.optimizations)
        else:
            self._fpoptimizer = None

        # mangled-name cache, keyed by original symbol, for this model
        self._mangled_names: dict[sp.Symbol, str] = {}
        # mangled names already assigned, for collision detection
        self._mangled_name_set: set[str] = set()

    def mangle_identifier(self, symbol: sp.Symbol) -> str:
        """Mangle the identifier for `symbol`, deduplicating against prior
        results for this model.

        The same symbol always yields the same output; distinct symbols
        never yield the same output -- keyed on the full symbol (name *and*
        assumptions), not just its name, since two AMICI-internal symbols
        can otherwise legitimately share a name with a differently-created
        (e.g. user-entity) symbol of the same name (#3240).
        """
        if (cached := self._mangled_names.get(symbol)) is not None:
            return cached
        base = mangled = _mangle(symbol.name)
        n = 2
        while mangled in self._mangled_name_set:
            mangled = f"{base}{n}"
            n += 1
        self._mangled_names[symbol] = mangled
        self._mangled_name_set.add(mangled)
        return mangled

    def _print_Symbol(self, expr: sp.Symbol) -> str:
        if expr == amici_time_symbol:
            return "t"
        return self.mangle_identifier(expr)

    def doprint(self, expr: sp.Expr, assign_to: str | None = None) -> str:
        if self._fpoptimizer:
            if isinstance(expr, list):
                expr = list(map(self._fpoptimizer, expr))
            else:
                expr = self._fpoptimizer(expr)

        try:
            # floating point
            code = super().doprint(expr, assign_to)
            code = re.sub(r"(^|\W)M_PI(\W|$)", r"\1amici::pi\2", code)

            return code
        except TypeError as e:
            raise ValueError(
                f'Encountered unsupported function in expression "{expr}"'
            ) from e

    def _print_min_max(self, expr, cpp_fun: str, sympy_fun):
        # C++ doesn't like mixing int and double for arguments for min/max,
        #  therefore, we just always convert to float
        args = [
            self._print(sp.Float(arg) if arg.is_number else arg)
            for arg in expr.args
        ]
        if len(expr.args) == 1:
            return args[0]
        if len(expr.args) == 2:
            return f"{self._ns}{cpp_fun}({args[0]}, {args[1]})"
        return f"{self._ns}{cpp_fun}({get_initializer_list(args)})"

    def _print_Min(self, expr):
        from sympy.functions.elementary.miscellaneous import Min

        return self._print_min_max(expr, "min", Min)

    def _print_Max(self, expr):
        from sympy.functions.elementary.miscellaneous import Max

        return self._print_min_max(expr, "max", Max)

    def _print_Infinity(self, expr):
        return "std::numeric_limits<double>::infinity()"

    def _print_NegativeInfinity(self, expr):
        return "-std::numeric_limits<double>::infinity()"

    def _print_ComplexInfinity(self, expr):
        # Since -zoo==+zoo, expressions containing ComplexInfinity may yield
        # unexpected results as compared to IEEE 754.

        warnings.warn(
            "Expression contains ComplexInfinity. "
            "This may point to a bug in the model and may result in "
            "incorrect simulation results. "
            "A possible cause is a division by zero in the model equations, "
            "which is supported by IEEE 754, but not by sympy."
            "Please check your model equations for potential issues.",
            stacklevel=1,
        )

        return "std::numeric_limits<double>::infinity()"

    def _get_sym_lines_array(
        self,
        equations: sp.Matrix,
        variable: str,
        indent_level: int,
        indices: Sequence[int] | None = None,
    ) -> list[str]:
        """
        Generate C++ code for assigning symbolic terms in symbols to C++ array
        `variable`.

        :param equations:
            vectors of symbolic expressions
        :param variable:
            name of the C++ array to assign to
        :param indent_level:
            indentation level (number of leading blanks)
        :param indices:
            List of custom indices corresponding to entries in `equations`.
            If `None`, the indices will be 0..(N-1).

        :return:
            C++ code as list of lines
        """
        if indices is None:
            indices = range(len(equations))

        if self.extract_cse:
            # placeholder names only, never seen by users -- must still be
            # valid identifiers since they go through the usual declare-and-
            # store printing, so no `[`/`]` as in the actual array access
            placeholder_symbols = sp.Matrix(
                [sp.Symbol(f"{variable}{index}") for index in indices]
            )
            res = self._get_output_declarations(
                symbols=placeholder_symbols,
                variable=variable,
                indent_level=indent_level,
                indices=indices,
            ) + self._get_sym_lines_symbols(
                symbols=placeholder_symbols,
                equations=equations,
                variable=variable,
                indent_level=indent_level,
                indices=indices,
            )
            # make compound statement so that extracted subexpressions are
            #  scoped locally and can be used in switch-cases
            indent = " " * indent_level
            return [
                f"{indent}{{",
                *(f"{indent}{l}" for l in res),
                f"{indent}}}",
            ]

        return [
            " " * indent_level + f"{variable}[{index}] = {self.doprint(math)};"
            for index, math in zip(indices, equations, strict=True)
            if math not in [0, 0.0]
        ]

    def _get_output_declarations(
        self,
        symbols: sp.Matrix,
        variable: str,
        indent_level: int,
        indices: Sequence[int] | None = None,
    ) -> list[str]:
        """
        Declare a named, non-`const` reference bound to each entry of
        `variable`, for `_get_sym_lines_symbols` to assign through.

        :param symbols:
            vectors of symbols that equations are assigned to

        :param variable:
            name of the C++ array to assign to

        :param indent_level:
            indentation level (number of leading blanks)

        :param indices:
            Optional custom indices corresponding to entries in `symbols`.

        :return:
            C++ code as list of lines
        """
        if indices is None:
            indices = range(len(symbols))
        else:
            assert len(indices) == len(symbols)

        indent = " " * indent_level
        return [
            f"{indent}realtype &{self.doprint(sym)} = {variable}[{index}];"
            for index, sym in zip(indices, symbols, strict=True)
        ]

    def _get_sym_lines_symbols(
        self,
        symbols: sp.Matrix,
        equations: sp.Matrix,
        variable: str,
        indent_level: int,
        indices: Sequence[int] | None = None,
    ) -> list[str]:
        """
        Generate C++ code assigning each entry's expression through the
        reference declared for it by `_get_output_declarations` (which
        must be called first, with the same arguments).

        :param symbols:
            vectors of symbols that equations are assigned to

        :param equations:
            vectors of expressions

        :param variable:
            name of the C++ array to assign to

        :param indent_level:
            indentation level (number of leading blanks)

        :param indices:
            Optional custom indices corresponding to entries in `symbols`.

        :return:
            C++ code as list of lines
        """
        assert len(symbols) == len(equations)
        if indices is None:
            indices = range(len(symbols))
        else:
            assert len(indices) == len(symbols)

        indent = " " * indent_level

        if self.extract_cse:
            # Extract common subexpressions
            cse_sym_prefix = "amici_cse_"
            symbol_generator = numbered_symbols(
                cls=sp.Symbol, prefix=cse_sym_prefix
            )
            replacements, reduced_exprs = sp.cse(
                equations,
                symbols=symbol_generator,
                order="none",
                list=False,
            )
            if replacements:
                # we need toposort to handle the dependencies of extracted
                #  subexpressions
                expr_dict = dict(
                    itertools.chain(
                        zip(symbols, reduced_exprs, strict=True), replacements
                    )
                )
                sorted_symbols = toposort(
                    {
                        identifier: {
                            s
                            for s in definition.free_symbols
                            if s in expr_dict
                        }
                        for (identifier, definition) in expr_dict.items()
                    }
                )
                symbol_to_idx = {
                    sym: idx for idx, sym in zip(indices, symbols, strict=True)
                }

                def format_line(symbol: sp.Symbol) -> str | None:
                    math = expr_dict[symbol]
                    if str(symbol).startswith(cse_sym_prefix):
                        return (
                            f"{indent}const realtype "
                            f"{self.doprint(symbol)} = "
                            f"{self.doprint(math)};"
                        )
                    if math in [0, 0.0]:
                        return None
                    math_str = self.doprint(math).replace("\n", "\n" + indent)
                    idx = symbol_to_idx[symbol]
                    return (
                        f"{indent}{self.doprint(symbol)} = {math_str};"
                        f"  // {variable}[{idx}]"
                    )

                return [
                    line
                    for symbol_group in sorted_symbols
                    for symbol in sorted(symbol_group, key=str)
                    if (line := format_line(symbol))
                ]

        lines = []
        for index, sym, math in zip(indices, symbols, equations, strict=True):
            if math in [0, 0.0]:
                continue
            math_str = self.doprint(math).replace("\n", "\n" + indent)
            lines.append(
                f"{indent}{self.doprint(sym)} = {math_str};"
                f"  // {variable}[{index}]"
            )
        return lines

    @staticmethod
    def print_bool(expr) -> str:
        """Print the boolean value of the given expression"""
        return "true" if bool(expr) else "false"


def get_switch_statement(
    condition: str,
    cases: dict[int, list[str]],
    indentation_level: int | None = 0,
    indentation_step: str | None = " " * 4,
) -> list[str]:
    """
    Generate code for a C++ switch statement.

    Generate code for a C++ switch statement with a ``break`` after each case.

    :param condition:
        Condition for switch

    :param cases:
        Cases as dict with expressions as keys and statement as
        list of strings

    :param indentation_level:
        indentation level

    :param indentation_step:
        indentation whitespace per level

    :return:
        Code for switch expression as list of strings
    """
    if not cases:
        return []

    indent0 = indentation_level * indentation_step
    indent1 = (indentation_level + 1) * indentation_step
    indent2 = (indentation_level + 2) * indentation_step

    # try to find redundant statements and collapse those cases
    # map statements to case expressions
    cases_map: dict[tuple[str, ...], list[str]] = {}
    for expression, statements in cases.items():
        if statements:
            statement_code = tuple(
                [
                    *(f"{indent2}{statement}" for statement in statements),
                    f"{indent2}break;",
                ]
            )
            case_code = f"{indent1}case {expression}:"

            cases_map[statement_code] = cases_map.get(statement_code, []) + [
                case_code
            ]

    if not cases_map:
        return []

    return [
        f"{indent0}switch({condition}) {{",
        *(
            code
            for codes in cases_map.items()
            for code in itertools.chain.from_iterable(reversed(codes))
        ),
        indent0 + "}",
    ]


def csc_matrix(
    matrix: sp.Matrix,
    rownames: list[sp.Symbol],
    colnames: list[sp.Symbol],
    identifier: int | None = 0,
    pattern_only: bool | None = False,
) -> tuple[list[int], list[int], sp.Matrix, list[sp.Symbol], sp.Matrix]:
    """
    Generates the sparse symbolic identifiers, symbolic identifiers,
    sparse matrix, column pointers and row values for a symbolic
    variable

    :param matrix:
        dense matrix to be sparsified

    :param rownames:
        ids of the variable of which the derivative is computed (assuming
        matrix is the jacobian)

    :param colnames:
        ids of the variable with respect to which the derivative is computed
        (assuming matrix is the jacobian)

    :param identifier:
        additional identifier that gets appended to symbol names to
        ensure their uniqueness in outer loops

    :param pattern_only:
        flag for computing sparsity pattern without whole matrix

    :return:
        symbol_col_ptrs, symbol_row_vals, sparse_list, symbol_list,
        sparse_matrix
    """
    idx = 0
    nrows, ncols = matrix.shape

    if not pattern_only:
        sparse_matrix = sp.zeros(nrows, ncols)
    symbol_list = []
    sparse_list = []
    symbol_col_ptrs = []
    symbol_row_vals = []

    for col in range(ncols):
        symbol_col_ptrs.append(idx)
        for row in range(nrows):
            if matrix[row, col].is_zero:
                continue

            symbol_row_vals.append(row)
            idx += 1
            symbol_name = f"d{rownames[row].name}_d{colnames[col].name}"
            if identifier:
                symbol_name += f"_{identifier}"
            symbol_list.append(symbol_with_assumptions(symbol_name))
            if pattern_only:
                continue

            sparse_matrix[row, col] = symbol_with_assumptions(symbol_name)
            sparse_list.append(matrix[row, col])

    if idx == 0:
        symbol_col_ptrs = []  # avoid bad memory access for empty matrices
    else:
        symbol_col_ptrs.append(idx)

    if pattern_only:
        sparse_matrix = None
    else:
        sparse_list = sp.Matrix(sparse_list)

    return (
        symbol_col_ptrs,
        symbol_row_vals,
        sparse_list,
        symbol_list,
        sparse_matrix,
    )


def get_initializer_list(values: Iterable) -> str:
    """Generate C++ initializer list for given values.

    :param values:
        Values to be included in the initializer list.
        They will be converted to strings, assuming :func:`str` will yield
        valid C++ expressions.
    """
    return f"{{{', '.join(map(str, values))}}}"
