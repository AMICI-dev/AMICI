import pytest
import sympy as sp
from amici.exporters.sundials.cxxcodeprinter import AmiciCxxCodePrinter
from amici.testing import skip_on_valgrind
from sympy.codegen.rewriting import optims_c99


@skip_on_valgrind
def test_optimizations():
    """Check that AmiciCxxCodePrinter handles optimizations correctly."""
    try:
        old_optim = AmiciCxxCodePrinter.optimizations
        AmiciCxxCodePrinter.optimizations = optims_c99
        cp = AmiciCxxCodePrinter()
        assert "expm1" in cp.doprint(sp.sympify("exp(x) - 1"))
    finally:
        AmiciCxxCodePrinter.optimizations = old_optim


@skip_on_valgrind
def test_print_infinity():
    """Check that AmiciCxxCodePrinter prints infinity correctly."""
    from sympy.core.numbers import ComplexInfinity, Infinity, NegativeInfinity

    cp = AmiciCxxCodePrinter()
    assert cp.doprint(Infinity()) == "std::numeric_limits<double>::infinity()"
    assert (
        cp.doprint(NegativeInfinity())
        == "-std::numeric_limits<double>::infinity()"
    )
    assert (
        cp.doprint(-NegativeInfinity())
        == "std::numeric_limits<double>::infinity()"
    )
    assert (
        cp.doprint(-Infinity()) == "-std::numeric_limits<double>::infinity()"
    )
    assert cp.doprint(sp.oo) == "std::numeric_limits<double>::infinity()"
    assert cp.doprint(-sp.oo) == "-std::numeric_limits<double>::infinity()"

    with pytest.warns(UserWarning, match="contains ComplexInfinity"):
        assert (
            cp.doprint(ComplexInfinity())
            == "std::numeric_limits<double>::infinity()"
        )
        assert (
            cp.doprint(-ComplexInfinity())
            == "std::numeric_limits<double>::infinity()"
        )
        assert cp.doprint(sp.zoo) == "std::numeric_limits<double>::infinity()"
        assert cp.doprint(-sp.zoo) == "std::numeric_limits<double>::infinity()"


@skip_on_valgrind
def test_min_max():
    """Check that AmiciCxxCodePrinter prints min() and max() correctly."""
    a, b, c = sp.symbols("a b c")
    cp = AmiciCxxCodePrinter()
    assert cp.doprint(sp.Min(a)) == "a_"
    assert cp.doprint(sp.Max(a)) == "a_"
    assert cp.doprint(sp.Min(a, b)) == "std::min(a_, b_)"
    assert cp.doprint(sp.Max(a, b)) == "std::max(a_, b_)"
    assert cp.doprint(sp.Min(a, b, c)) == "std::min({a_, b_, c_})"
    assert cp.doprint(sp.Max(a, b, c)) == "std::max({a_, b_, c_})"


@skip_on_valgrind
def test_float_arithmetic():
    """
    Check that AmiciCxxCodePrinter produces code that uses float arithmetic.
    """
    cp = AmiciCxxCodePrinter()
    assert cp.doprint(sp.Rational(1, 2)) == "1.0/2.0"
    assert cp.doprint(sp.Integer(1) / sp.Integer(2)) == "1.0/2.0"


@skip_on_valgrind
def test_mangle_identifier():
    """Check identifier mangling for reserved/keyword-colliding names."""
    cp = AmiciCxxCodePrinter()

    # ordinary names just get a trailing underscore
    assert cp.mangle_identifier(sp.Symbol("STAT")) == "STAT_"

    # names already ending in `_` get a bare-letter marker, not `_v` --
    # regression test: `f"{name}_v"` would produce "k__v", reintroducing
    # exactly the `__` this is supposed to avoid
    assert cp.mangle_identifier(sp.Symbol("k_")) == "k_v"
    assert "__" not in cp.mangle_identifier(sp.Symbol("k_"))
    assert "__" not in cp.mangle_identifier(sp.Symbol("k__"))

    # an internal `__` run collapses to a single `_` before suffixing
    assert cp.mangle_identifier(sp.Symbol("my__species")) == "my_species_"

    # C++ keywords and stdlib macros mangle to something distinct from the
    # original token
    for name in (
        "int",
        "class",
        "template",
        "new",
        "for",
        "NULL",
        "EOF",
        "INFINITY",
        "NAN",
    ):
        mangled = cp.mangle_identifier(sp.Symbol(name))
        assert mangled != name
        assert "__" not in mangled

    # same input -> same output (cache hit, not recomputed)
    assert cp.mangle_identifier(sp.Symbol("STAT")) == cp.mangle_identifier(
        sp.Symbol("STAT")
    )

    # distinct inputs that collapse to the same base still get distinct,
    # `__`-free results
    r1 = cp.mangle_identifier(sp.Symbol("a__b"))
    r2 = cp.mangle_identifier(sp.Symbol("a_b"))
    assert r1 != r2
    assert "__" not in r1
    assert "__" not in r2


@skip_on_valgrind
def test_mangle_identifier_distinguishes_symbols_with_same_name():
    """Two distinct symbols that merely happen to share a `.name` (e.g. one
    is a plain user-entity symbol, the other one of AMICI's own
    assumption-free internal placeholders for e.g. a state derivative or a
    spline) must never be mangled to the same C++ identifier -- otherwise
    two logically different quantities silently collide into a single
    declared local (#3240)."""
    cp = AmiciCxxCodePrinter()

    same_name_plain = sp.Symbol("dpdt")
    same_name_real = sp.Symbol("dpdt", real=True)
    assert same_name_plain != same_name_real  # sanity check on the premise

    r1 = cp.mangle_identifier(same_name_plain)
    r2 = cp.mangle_identifier(same_name_real)
    assert r1 != r2

    # each symbol is still cached by its own identity: repeated mangling of
    # the *same* object doesn't drift or get reassigned to the other's slot
    assert cp.mangle_identifier(same_name_plain) == r1
    assert cp.mangle_identifier(same_name_real) == r2


@skip_on_valgrind
def test_extract_cse():
    """Test extraction of common subexpressions."""
    cp = AmiciCxxCodePrinter()
    cp_cse = AmiciCxxCodePrinter(extract_cse=True)

    a, b, c = sp.symbols("a b c")
    x1, x2, x3 = sp.symbols("x1 x2 x3")

    syms = sp.Matrix([x1, x2, x3])
    eqs = sp.Matrix([a * b * c, a * b, a * b * c + a])

    # every output gets a reference bound to its array slot, declared
    # once, so the assignment itself can use the readable name as its LHS
    expected_decl = [
        "  realtype &x1_ = x[0];",
        "  realtype &x2_ = x[1];",
        "  realtype &x3_ = x[2];",
    ]
    expected = [
        "  x1_ = a_*b_*c_;  // x[0]",
        "  x2_ = a_*b_;  // x[1]",
        "  x3_ = a_*b_*c_ + a_;  // x[2]",
    ]
    expected_cse = [
        "  const realtype amici_cse_0_ = a_*b_;",
        "  const realtype amici_cse_1_ = amici_cse_0_*c_;",
        "  x2_ = amici_cse_0_;  // x[1]",
        "  x1_ = amici_cse_1_;  // x[0]",
        "  x3_ = a_ + amici_cse_1_;  // x[2]",
    ]

    assert expected_decl == cp._get_output_declarations(
        symbols=syms, variable="x", indent_level=2
    )
    assert expected == cp._get_sym_lines_symbols(
        symbols=syms, equations=eqs, variable="x", indent_level=2
    )
    assert expected_cse == cp_cse._get_sym_lines_symbols(
        symbols=syms, equations=eqs, variable="x", indent_level=2
    )

    expected = [
        "  x[0] = a_*b_*c_;",
        "  x[1] = a_*b_;",
        "  x[2] = a_*b_*c_ + a_;",
    ]
    expected_cse = [
        "  {",
        "    realtype &x0_ = x[0];",
        "    realtype &x1_ = x[1];",
        "    realtype &x2_ = x[2];",
        "    const realtype amici_cse_0_ = a_*b_;",
        "    const realtype amici_cse_1_ = amici_cse_0_*c_;",
        "    x1_ = amici_cse_0_;  // x[1]",
        "    x0_ = amici_cse_1_;  // x[0]",
        "    x2_ = a_ + amici_cse_1_;  // x[2]",
        "  }",
    ]
    assert expected == cp._get_sym_lines_array(
        equations=eqs, variable="x", indent_level=2
    )
    assert expected_cse == cp_cse._get_sym_lines_array(
        equations=eqs, variable="x", indent_level=2
    )


@skip_on_valgrind
def test_sym_lines_symbols_output_reference_reused():
    """A later equation can read an earlier output through its declared
    reference (e.g. `w`'s dynamic expressions depending on an earlier `w`
    entry) -- and, since the reference doesn't depend on any computed
    value, this works even across control-flow boundaries the equations
    themselves don't share (see the `include_static`/dynamic split in
    `DEExporter._get_function_body`)."""
    cp = AmiciCxxCodePrinter()
    a, b = sp.symbols("a b")
    y1, y2 = sp.symbols("y1 y2")
    syms = sp.Matrix([y1, y2])
    eqs = sp.Matrix([a * b, y1 + 1])  # y2's equation references y1

    expected_decl = [
        "  realtype &y1_ = x[0];",
        "  realtype &y2_ = x[1];",
    ]
    expected = [
        "  y1_ = a_*b_;  // x[0]",
        "  y2_ = y1_ + 1;  // x[1]",
    ]
    assert expected_decl == cp._get_output_declarations(
        symbols=syms, variable="x", indent_level=2
    )
    assert expected == cp._get_sym_lines_symbols(
        symbols=syms, equations=eqs, variable="x", indent_level=2
    )
