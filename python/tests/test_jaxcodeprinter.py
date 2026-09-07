"""Tests for the JAX code printer, in particular identifier mangling."""

import keyword

import pytest
import sympy as sp
from amici.testing import skip_on_valgrind

pytest.importorskip("jax")

from amici.exporters.jax.jaxcodeprinter import (  # noqa: E402
    AmiciJaxCodePrinter,
    generic_measurement_symbol,
)
from amici.importers.utils import (  # noqa: E402
    amici_time_symbol,
    symbol_with_assumptions,
)


@skip_on_valgrind
def test_mangle_identifier():
    """Every model entity id printed into the generated module is mangled,
    so it can't collide with anything that module already uses."""
    cp = AmiciJaxCodePrinter()

    # ordinary names just get a trailing underscore
    assert cp.doprint(sp.Symbol("STAT")) == "STAT_"

    # names already ending in `_` get a bare-letter marker, and internal
    # `__` runs collapse -- as for the C++ backend, so the two backends
    # name the same entity the same way
    assert cp.doprint(sp.Symbol("k_")) == "k_v"
    assert cp.doprint(sp.Symbol("my__species")) == "my_species_"

    # nothing the generated module itself uses can be produced: amici's
    # fixed array-parameter names, the generated methods' other arguments
    # and locals, the module-level imports, and Python keywords all mangle
    # to something else
    for name in (
        *("x", "p", "k", "h", "w", "y"),
        *("tcl", "op", "np", "my", "iy", "args", "self"),
        *("jnp", "jr", "eqx", "oo", "safe_log", "safe_div", "JAXModel"),
        *keyword.kwlist,
    ):
        # created the way the importers create model entities
        mangled = cp.doprint(symbol_with_assumptions(name))
        assert mangled != name
        assert not keyword.iskeyword(mangled)
        assert mangled.isidentifier()

    # distinct inputs that collapse to the same base still get distinct
    # results
    assert cp.doprint(sp.Symbol("a__b")) != cp.doprint(sp.Symbol("a_b"))

    # same input -> same output
    assert cp.doprint(sp.Symbol("STAT")) == "STAT_"


@skip_on_valgrind
def test_fixed_symbols_are_not_mangled():
    """The two symbols that denote a fixed argument of the generated
    functions rather than a model entity keep the name that argument has --
    while a model entity that merely shares their name doesn't."""
    cp = AmiciJaxCodePrinter()

    assert cp.doprint(amici_time_symbol) == "t"
    assert cp.doprint(generic_measurement_symbol) == "my"

    # a model entity named `my` is a different symbol and gets a mangled
    # name of its own (an entity named `t` can't occur -- `t` is the one
    # name still renamed at model level, see RESERVED_SYMBOLS)
    entity_my = sp.Symbol("my", real=True)
    assert entity_my != generic_measurement_symbol  # premise
    assert cp.doprint(entity_my) == "my_"


@skip_on_valgrind
def test_assignment_targets_are_mangled_consistently():
    """An assignment's left-hand side has to use the same identifier the
    symbol gets inside an expression -- otherwise the generated code reads
    from the unmangled (and possibly shadowed, or syntactically invalid)
    name."""
    cp = AmiciJaxCodePrinter()

    # `class` is a Python keyword, `jnp` one of the module's imports
    class_, jnp_ = sp.Symbol("class"), sp.Symbol("jnp")
    lines = cp._get_sym_lines(
        sp.Matrix([class_, jnp_]),
        sp.Matrix([sp.Float(2.0), 3 * class_]),
        indent_level=0,
    )

    assert lines == [
        f"{cp.doprint(class_)} = 2.0",
        f"{cp.doprint(jnp_)} = 3*{cp.doprint(class_)}",
    ]
    compile("\n".join(lines), "<generated>", "exec")
