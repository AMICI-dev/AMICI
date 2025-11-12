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
