import sympy as sp
from amici.cxxcodeprinter import AmiciCxxCodePrinter
from sympy.codegen.rewriting import optims_c99


def test_optimizations():
    """Check that AmiciCxxCodePrinter handles optimizations correctly."""
    try:
        old_optim = AmiciCxxCodePrinter.optimizations
        AmiciCxxCodePrinter.optimizations = optims_c99
        cp = AmiciCxxCodePrinter()
        assert "expm1" in cp.doprint(sp.sympify("exp(x) - 1"))
    finally:
        AmiciCxxCodePrinter.optimizations = old_optim
