"""Jax code generation"""

import itertools
import math
import re
from collections.abc import Iterable
from logging import warning

import sympy as sp
from sympy.core.function import UndefinedFunction
from sympy.printing.numpy import NumPyPrinter
from toposort import toposort

from amici.exporters._mangling import IdentifierMangler
from amici.importers.utils import amici_time_symbol

#: The generic measurement symbol that observable-specific measurement
#: symbols are substituted by before code generation. Unlike every other
#: symbol printed here, it denotes a fixed argument of the generated
#: ``_nllh`` (see ``jax.template.py``), not a model entity.
generic_measurement_symbol = sp.Symbol("my")


class AmiciJaxCodePrinter(NumPyPrinter):
    """JAX code printer"""

    #: Symbols denoting one of the fixed arguments of the generated
    #: functions rather than a model entity, mapped to the argument name
    #: they have to be printed as (see ``jax.template.py``). Keyed on the
    #: full symbol, so a model entity that merely shares the name is
    #: unaffected (and gets a mangled name of its own, as any other).
    _fixed_symbol_names = {
        amici_time_symbol: "t",
        generic_measurement_symbol: "my",
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # identifier mangling for this model, shared with the C++ backend
        self._mangler = IdentifierMangler()

    def mangle_identifier(self, symbol: sp.Symbol) -> str:
        """Mangle the identifier for `symbol` into a safe, unique Python
        local variable name.

        Every model entity id printed into the generated module -- as the
        target of an assignment as much as inside an expression -- has to
        go through here. The generated functions destructure their array
        arguments into per-entry locals (``x_, y_ = x``), so an unmangled
        entity id would shadow the argument it was unpacked from (or the
        module's imports, or a Python keyword), silently corrupting every
        later use of it in that function.

        See :meth:`amici.exporters._mangling.IdentifierMangler.mangle`.
        """
        return self._mangler.mangle(symbol)

    def _print_Symbol(self, expr: sp.Symbol) -> str:
        if (fixed_name := self._fixed_symbol_names.get(expr)) is not None:
            return fixed_name
        return self.mangle_identifier(expr)

    def _print_Float(self, expr: sp.Float) -> str:
        # sympy's string printer emits 15 significant digits, which is not
        # enough to round-trip a double. The lost bits change results wherever
        # a value is compared exactly, e.g. a parameter with the value `pi`
        # would no longer test equal to `pi` (SBML test suite case 00958).
        # Emit the shortest representation that recovers the exact double
        # instead (the SUNDIALS backend prints 17 digits for the same reason).
        value = float(expr)
        if not math.isfinite(value):
            # outside the double range; `repr` would emit a bare `inf`, which
            # is not a valid literal in the generated module
            return super()._print_Float(expr)
        return repr(value)

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
        if not isinstance(expr.func, UndefinedFunction):
            return super()._print_Function(expr)

        nn_name = expr.func.__name__
        input_args = expr.args[:-1]
        output_idx = expr.args[-1]
        forward_arg = self._build_nn_forward_arg(input_args)
        return f"self.nns['{nn_name}'].forward({forward_arg})[{output_idx}]"

    def _build_nn_forward_arg(self, input_args) -> str:
        array_prefix = "_nn_array_"
        if not any(
            hasattr(a, "name") and a.name.startswith(array_prefix)
            for a in input_args
        ):
            return f"jnp.array([{', '.join(self.doprint(a) for a in input_args)}])"

        groups = self._group_nn_inputs(input_args, array_prefix)
        return groups[0] if len(groups) == 1 else f"[{', '.join(groups)}]"

    def _group_nn_inputs(self, input_args, array_prefix: str) -> list[str]:
        groups = []
        scalar_group = []
        for a in input_args:
            if hasattr(a, "name") and a.name.startswith(array_prefix):
                if scalar_group:
                    groups.append(f"jnp.array([{', '.join(scalar_group)}])")
                    scalar_group = []
                petab_id = a.name[len(array_prefix) :]
                groups.append(
                    f"self._array_inputs['{petab_id}'][self._array_input_index]"
                )
            else:
                scalar_group.append(self.doprint(a))
        if scalar_group:
            groups.append(f"jnp.array([{', '.join(scalar_group)}])")
        return groups

    def _print_Max(self, expr: sp.Expr) -> str:
        """
        Print the max function, replacing it with jnp.max.
        """
        return f"jnp.max({_jnp_array_str(expr.args, self)})"

    def _print_Min(self, expr: sp.Expr) -> str:
        """
        Print the min function, replacing it with jnp.min.
        """
        return f"jnp.min({_jnp_array_str(expr.args, self)})"

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
        symbols = list(symbols)
        equations = list(equations)

        # Extract common subexpressions before emitting the equations. Some
        # generated blocks (e.g. analytic steady-state initial values) repeat
        # large subterms many times; JAX's autodiff differentiates every copy,
        # producing a huge reverse-mode graph whose XLA compilation dominates
        # runtime. Factoring shared terms into temporaries keeps the
        # differentiated graph compact (and mirrors the SUNDIALS backend).
        replacements, reduced = sp.cse(
            equations,
            symbols=sp.numbered_symbols(prefix="_amici_cse_", cls=sp.Symbol),
            order="none",
            list=False,
        )
        if not replacements:
            return [
                f"{indent}{self._print_assignment_target(s)} = "
                f"{self.doprint(e)}"
                for s, e in zip(symbols, reduced)
            ]

        # Assignment blocks are sequential: a target symbol (e.g. a `w` entry)
        # may be referenced by later targets *and* by extracted temporaries, so
        # temporaries cannot simply be emitted first. Toposort the combined set
        # of targets and temporaries by their mutual dependencies so every
        # symbol is defined before use (mirrors the SUNDIALS code printer).
        expr_dict = dict(
            itertools.chain(zip(symbols, reduced, strict=True), replacements)
        )
        dependencies = {
            identifier: {s for s in definition.free_symbols if s in expr_dict}
            for identifier, definition in expr_dict.items()
        }
        return [
            f"{indent}{self._print_assignment_target(sym)} = "
            f"{self.doprint(expr_dict[sym])}"
            for group in toposort(dependencies)
            for sym in sorted(group, key=str)
        ]

    def _print_assignment_target(self, symbol: sp.Symbol | str) -> str:
        """Print `symbol` as the left-hand side of an assignment.

        Same identifier as the symbol gets inside an expression -- the two
        sides have to agree, which they only do if both go through the
        printer.
        """
        return (
            self.doprint(symbol)
            if isinstance(symbol, sp.Basic)
            else str(symbol)
        )


def _jnp_array_str(
    array, code_printer: AmiciJaxCodePrinter | None = None
) -> str:
    printer = code_printer or AmiciJaxCodePrinter()
    elems = ", ".join(
        printer.doprint(s) if isinstance(s, sp.Basic) else str(s)
        for s in array
    )

    return f"jnp.array([{elems}])"
