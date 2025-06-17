import sympy as sp
from amici.cxxcodeprinter import AmiciCxxCodePrinter
from sympy.codegen.rewriting import optims_c99
from amici.testing import skip_on_valgrind


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
    from sympy.core.numbers import NegativeInfinity, Infinity, ComplexInfinity

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
