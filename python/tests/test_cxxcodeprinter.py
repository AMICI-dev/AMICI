from amici.cxxcodeprinter import AmiciCxxCodePrinter
from sympy.codegen.rewriting import optims_c99
import sympy as sp


def test_optimizations():
    """Check that AmiciCxxCodePrinter handles optimizations correctly."""
    try:
        old_optim = AmiciCxxCodePrinter.optimizations
        AmiciCxxCodePrinter.optimizations = optims_c99
        cp = AmiciCxxCodePrinter()
        assert "expm1" in cp.doprint(sp.sympify("exp(x) - 1"))
    except Exception:
        raise
    finally:
        AmiciCxxCodePrinter.optimizations = old_optim
