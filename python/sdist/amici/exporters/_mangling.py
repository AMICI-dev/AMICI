"""Turning model entity ids into identifiers that are safe in generated code.

Both backends emit model entity ids as identifiers in the code they
generate. A model entity id is only constrained by the source format
(any SBML ``SId``, any PySB component name), so it may well collide with
something the target language or the generated code itself already uses:
a keyword, a standard-library macro, one of the fixed argument names the
generated functions are written with, or another entity id that isn't
distinguishable after mangling.

Passing every printed identifier through a single choke point per model
(:class:`IdentifierMangler`) guarantees none of that can happen, without
the model itself having to know anything about the target language.
"""

import re

import sympy as sp

__all__ = ["IdentifierMangler", "mangle_name"]


def mangle_name(name: str) -> str:
    """Make a model-derived identifier safe as a local variable name in
    generated code.

    Appends `_` so it can never equal a real keyword/macro (in either
    target language: no C++ keyword, C standard-library macro, Python
    keyword or builtin ends in an underscore, and neither does any of the
    fixed argument/variable names of the generated functions). Collapses
    any SBML-legal `__` run first, and uses `v` instead of `_` as the
    marker for names already ending in `_`, so the result never contains
    `__` either (reserved in C++; class-private in Python).

    This is not injective and does not guarantee a unique result on its
    own (e.g. `"a__b"` and `"a_b"` both collapse to the same string) --
    callers that need uniqueness across a whole model handle that
    separately (see :class:`IdentifierMangler`).
    """
    name = re.sub(r"_{2,}", "_", name)
    return f"{name}v" if name.endswith("_") else f"{name}_"


class IdentifierMangler:
    """Assigns each symbol of one model a unique, safe identifier for the
    generated code.

    One instance per generated model: the mangled name of a given symbol
    has to be the same everywhere it occurs, and must not be handed to any
    other symbol of that model.
    """

    def __init__(self):
        # mangled-name cache, keyed by original symbol, for this model
        self._mangled_names: dict[sp.Symbol, str] = {}
        # mangled names already assigned, for collision detection
        self._mangled_name_set: set[str] = set()

    def mangle(self, symbol: sp.Symbol) -> str:
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
        base = mangled = mangle_name(symbol.name)
        n = 2
        while mangled in self._mangled_name_set:
            mangled = f"{base}{n}"
            n += 1
        self._mangled_names[symbol] = mangled
        self._mangled_name_set.add(mangled)
        return mangled
