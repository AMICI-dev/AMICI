"""C++ code generation"""
import re
from typing import List, Optional, Tuple, Dict

import sympy as sp
from sympy.printing.cxx import CXX11CodePrinter


class AmiciCxxCodePrinter(CXX11CodePrinter):
    """C++ code printer"""

    def __init__(self):
        super().__init__()

    def doprint(self, expr: sp.Expr, assign_to: Optional[str] = None) -> str:
        try:
            code = super().doprint(expr, assign_to)
            code = re.sub(r'(^|\W)M_PI(\W|$)', r'\1amici::pi\2', code)

            return code
        except TypeError as e:
            raise ValueError(
                f'Encountered unsupported function in expression "{expr}"'
            ) from e

    def _print_min_max(self, expr, cpp_fun: str, sympy_fun):
        # C++ doesn't like mixing int and double for arguments for min/max,
        #  therefore, we just always convert to float
        arg0 = sp.Float(expr.args[0]) if expr.args[0].is_number \
            else expr.args[0]
        if len(expr.args) == 1:
            return self._print(arg0)
        return "%s%s(%s, %s)" % (self._ns, cpp_fun, self._print(arg0),
                                 self._print(sympy_fun(*expr.args[1:])))

    def _print_Min(self, expr):
        from sympy.functions.elementary.miscellaneous import Min
        return self._print_min_max(expr, "min", Min)

    def _print_Max(self, expr):
        from sympy.functions.elementary.miscellaneous import Max
        return self._print_min_max(expr, "max", Max)

    def _get_sym_lines_array(
            self,
            equations: sp.Matrix,
            variable: str,
            indent_level: int
    ) -> List[str]:
        """
        Generate C++ code for assigning symbolic terms in symbols to C++ array
        `variable`.

        :param equations:
            vectors of symbolic expressions

        :param variable:
            name of the C++ array to assign to

        :param indent_level:
            indentation level (number of leading blanks)

        :return:
            C++ code as list of lines
        """
        return [
            ' ' * indent_level + f'{variable}[{index}] = '
                                 f'{self.doprint(math)};'
            for index, math in enumerate(equations)
            if math not in [0, 0.0]
        ]

    def _get_sym_lines_symbols(
            self, symbols: sp.Matrix,
            equations: sp.Matrix,
            variable: str,
            indent_level: int
    ) -> List[str]:
        """
        Generate C++ code for where array elements are directly replaced with
        their corresponding macro symbol

        :param symbols:
            vectors of symbols that equations are assigned to

        :param equations:
            vectors of expressions

        :param variable:
            name of the C++ array to assign to, only used in comments

        :param indent_level:
            indentation level (number of leading blanks)

        :return:
            C++ code as list of lines
        """
        return [
            f'{" " * indent_level}{sym} = {self.doprint(math)};'
            f'  // {variable}[{index}]'.replace('\n',
                                                '\n' + ' ' * indent_level)
            for index, (sym, math) in enumerate(zip(symbols, equations))
            if math not in [0, 0.0]
        ]

    def csc_matrix(
            self,
            matrix: sp.Matrix,
            rownames: List[sp.Symbol],
            colnames: List[sp.Symbol],
            identifier: Optional[int] = 0,
            pattern_only: Optional[bool] = False
    ) -> Tuple[
        List[int], List[int], sp.Matrix, List[str], sp.Matrix
    ]:
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
                if matrix[row, col] == 0:
                    continue

                symbol_row_vals.append(row)
                idx += 1
                symbol_name = f'd{self.doprint(rownames[row])}' \
                              f'_d{self.doprint(colnames[col])}'
                if identifier:
                    symbol_name += f'_{identifier}'
                symbol_list.append(symbol_name)
                if pattern_only:
                    continue

                sparse_matrix[row, col] = sp.Symbol(symbol_name, real=True)
                sparse_list.append(matrix[row, col])

        if idx == 0:
            symbol_col_ptrs = []  # avoid bad memory access for empty matrices
        else:
            symbol_col_ptrs.append(idx)

        if pattern_only:
            sparse_matrix = None
        else:
            sparse_list = sp.Matrix(sparse_list)

        return symbol_col_ptrs, symbol_row_vals, sparse_list, symbol_list, \
               sparse_matrix


def get_switch_statement(condition: str, cases: Dict[int, List[str]],
                         indentation_level: Optional[int] = 0,
                         indentation_step: Optional[str] = ' ' * 4):
    """
    Generate code for switch statement

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
    lines = []

    if not cases:
        return lines

    for expression, statements in cases.items():
        if statements:
            lines.append((indentation_level + 1) * indentation_step
                         + f'case {expression}:')
            for statement in statements:
                lines.append((indentation_level + 2) * indentation_step
                             + statement)
            lines.append((indentation_level + 2) * indentation_step + 'break;')

    if lines:
        lines.insert(0, indentation_level * indentation_step
                     + f'switch({condition}) {{')
        lines.append(indentation_level * indentation_step + '}')

    return lines

