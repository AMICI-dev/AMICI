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
    assert cp.doprint(sp.Min(a)) == "a"
    assert cp.doprint(sp.Max(a)) == "a"
    assert cp.doprint(sp.Min(a, b)) == "std::min(a, b)"
    assert cp.doprint(sp.Max(a, b)) == "std::max(a, b)"
    assert cp.doprint(sp.Min(a, b, c)) == "std::min({a, b, c})"
    assert cp.doprint(sp.Max(a, b, c)) == "std::max({a, b, c})"


@skip_on_valgrind
def test_float_arithmetic():
    """
    Check that AmiciCxxCodePrinter produces code that uses float arithmetic.
    """
    cp = AmiciCxxCodePrinter()
    assert cp.doprint(sp.Rational(1, 2)) == "1.0/2.0"
    assert cp.doprint(sp.Integer(1) / sp.Integer(2)) == "1.0/2.0"


@skip_on_valgrind
def test_extract_cse():
    """Test extraction of common subexpressions."""
    cp = AmiciCxxCodePrinter()
    cp_cse = AmiciCxxCodePrinter(extract_cse=True)

    a, b, c = sp.symbols("a b c")
    x1, x2, x3 = sp.symbols("x1 x2 x3")

    syms = sp.Matrix([x1, x2, x3])
    eqs = sp.Matrix([a * b * c, a * b, a * b * c + a])

    expected = [
        "  x1 = a*b*c;  // x[0]",
        "  x2 = a*b;  // x[1]",
        "  x3 = a*b*c + a;  // x[2]",
    ]

    expected_cse = [
        "  const realtype __amici_cse_0 = a*b;",
        "  const realtype __amici_cse_1 = __amici_cse_0*c;",
        "  x2 = __amici_cse_0;  // x[1]",
        "  x1 = __amici_cse_1;  // x[0]",
        "  x3 = __amici_cse_1 + a;  // x[2]",
    ]

    assert expected == cp._get_sym_lines_symbols(
        symbols=syms, equations=eqs, variable="x", indent_level=2
    )
    assert expected_cse == cp_cse._get_sym_lines_symbols(
        symbols=syms, equations=eqs, variable="x", indent_level=2
    )

    expected = ["  x[0] = a*b*c;", "  x[1] = a*b;", "  x[2] = a*b*c + a;"]
    expected_cse = [
        "  {",
        "    const realtype __amici_cse_0 = a*b;",
        "    const realtype __amici_cse_1 = __amici_cse_0*c;",
        "    x[1] = __amici_cse_0;  // x[1]",
        "    x[0] = __amici_cse_1;  // x[0]",
        "    x[2] = __amici_cse_1 + a;  // x[2]",
        "  }",
    ]
    assert expected == cp._get_sym_lines_array(
        equations=eqs, variable="x", indent_level=2
    )
    assert expected_cse == cp_cse._get_sym_lines_array(
        equations=eqs, variable="x", indent_level=2
    )
