"""Jax code generation"""

import re
from collections.abc import Iterable
from logging import warning

import sympy as sp
from sympy.core.function import UndefinedFunction
from sympy.printing.numpy import NumPyPrinter


def _jnp_array_str(array) -> str:
    elems = ", ".join(str(s) for s in array)

    return f"jnp.array([{elems}])"


class AmiciJaxCodePrinter(NumPyPrinter):
    """JAX code printer"""

    def doprint(self, expr: sp.Expr, assign_to: str | None = None) -> str:
        try:
            code = super().doprint(expr, assign_to)
            code = re.sub(r"numpy\.", r"jnp.", code)

            return code
        except TypeError as e:
            raise ValueError(
                f'Encountered unsupported function in expression "{expr}"'
            ) from e

    def _print_AmiciSpline(self, expr: sp.Expr) -> str:
        warning("Spline interpolation is support in JAX is untested")
        # FIXME: untested, where are spline nodes coming from anyways?
        return f'interp1d(time, {self.doprint(expr.args[2:])}, kind="cubic")'

    def _print_log(self, expr: sp.Expr) -> str:
        return f"safe_log({self.doprint(expr.args[0])})"

    def _print_Mul(self, expr: sp.Expr) -> str:
        numer, denom = expr.as_numer_denom()
        if denom == 1:
            return super()._print_Mul(expr)
        return f"safe_div({self.doprint(numer)}, {self.doprint(denom)})"

    def _print_Function(self, expr: sp.Expr) -> str:
        if isinstance(expr.func, UndefinedFunction):
            nn_name = expr.func.__name__
            input_args = expr.args[:-1]
            output_idx = expr.args[-1]

            # Check if any input args are array inputs (convention: _nn_array_<petab_id>)
            array_prefix = "_nn_array_"
            has_array = any(
                hasattr(a, "name") and a.name.startswith(array_prefix)
                for a in input_args
            )

            if has_array:
                # Build grouped input: scalars go in jnp.array([...]),
                # array symbols become self._array_inputs['<petab_id>']
                groups = []
                scalar_group = []
                for a in input_args:
                    if hasattr(a, "name") and a.name.startswith(array_prefix):
                        # Flush any accumulated scalars as a group
                        if scalar_group:
                            groups.append(
                                f"jnp.array([{', '.join(scalar_group)}])"
                            )
                            scalar_group = []
                        petab_id = a.name[len(array_prefix) :]
                        groups.append(
                            f"self._array_inputs['{petab_id}'][self._array_input_index]"
                        )
                    else:
                        scalar_group.append(self.doprint(a))
                if scalar_group:
                    groups.append(
                        f"jnp.array([{', '.join(scalar_group)}])"
                    )

                if len(groups) == 1:
                    forward_arg = groups[0]
                else:
                    forward_arg = f"[{', '.join(groups)}]"

                return (
                    f"self.nns['{nn_name}'].forward("
                    f"{forward_arg})[{output_idx}]"
                )

            return (
                f"self.nns['{nn_name}'].forward("
                f"jnp.array([{', '.join(self.doprint(a) for a in input_args)}])"
                f")[{output_idx}]"
            )
        else:
            return super()._print_Function(expr)

    def _print_Max(self, expr: sp.Expr) -> str:
        """
        Print the max function, replacing it with jnp.max.
        """
        return f"jnp.max({_jnp_array_str(expr.args)})"

    def _print_Min(self, expr: sp.Expr) -> str:
        """
        Print the min function, replacing it with jnp.min.
        """
        return f"jnp.min({_jnp_array_str(expr.args)})"

    def _get_sym_lines(
        self,
        symbols: sp.Matrix | Iterable[str],
        equations: sp.Matrix | Iterable[sp.Expr],
        indent_level: int,
    ) -> list[str]:
        """
        Generate C++ code for assigning symbolic terms in symbols to C++ array
        `variable`.

        :param equations:
            vectors of symbolic expressions

        :param symbols:
            names of the symbols to assign to

        :param indent_level:
            indentation level (number of leading blanks)

        :return:
            C++ code as list of lines
        """
        indent = " " * indent_level
        return [
            f"{indent}{s} = {self.doprint(e)}"
            for s, e in zip(symbols, equations)
        ]
