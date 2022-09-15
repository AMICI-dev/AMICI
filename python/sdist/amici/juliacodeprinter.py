"""Julia code generation"""
import itertools
import os
import re
from typing import List, Optional, Union, Iterable

import sympy as sp
from sympy.printing.julia import JuliaCodePrinter


class AmiciJuliaCodePrinter(JuliaCodePrinter):
    """Julia code printer"""

    def doprint(self, expr: sp.Expr, assign_to: Optional[str] = None) -> str:
        try:
            code = super().doprint(expr, assign_to)

            return code
        except TypeError as e:
            raise ValueError(
                f'Encountered unsupported function in expression "{expr}"'
            ) from e

    def _get_sym_lines(
            self,
            symbols: Union[Iterable[str], sp.Matrix],
            equations: sp.Matrix,
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
            Julia code as list of lines
        """
        indent = ' ' * indent_level
        return [
            f'{indent}{s} = {self.doprint(e)}'
            for s, e in zip(symbols, equations)
        ]
