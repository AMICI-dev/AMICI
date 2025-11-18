"""Tests related to the sympy_utils module."""

import sympy as sp
from amici._symbolic.sympy_utils import (
    _custom_pow_eval_derivative,
    _monkeypatched,
    _piecewise_to_minmax,
)
from amici.testing import skip_on_valgrind


@skip_on_valgrind
def test_monkeypatch():
    t = sp.Symbol("t")
    n = sp.Symbol("n")
    vals = [(t, 0), (n, 1)]

    # check that the removable singularity still exists
    assert (t**n).diff(t).subs(vals) is sp.nan

    # check that we can monkeypatch it out
    with _monkeypatched(
        sp.Pow, "_eval_derivative", _custom_pow_eval_derivative
    ):
        assert (t**n).diff(t).subs(vals) is not sp.nan

    # check that the monkeypatch is transient
    assert (t**n).diff(t).subs(vals) is sp.nan


@skip_on_valgrind
def test_rewrite_piecewise_minmax():
    """Test rewriting of piecewise min/max to sympy Min/Max functions."""
    x, y, z = sp.symbols("x y z")

    assert sp.Piecewise((x, x < y), (y, True)).replace(
        sp.Piecewise, _piecewise_to_minmax
    ) == sp.Min(x, y)
    assert sp.Piecewise((x, x <= y), (y, True)).replace(
        sp.Piecewise, _piecewise_to_minmax
    ) == sp.Min(x, y)
    assert sp.Piecewise((x, x > y), (y, True)).replace(
        sp.Piecewise, _piecewise_to_minmax
    ) == sp.Max(x, y)
    assert sp.Piecewise((x, x >= y), (y, True)).replace(
        sp.Piecewise, _piecewise_to_minmax
    ) == sp.Max(x, y)
    assert sp.Piecewise((x, y > x), (y, True)).replace(
        sp.Piecewise, _piecewise_to_minmax
    ) == sp.Min(x, y)
    assert sp.Piecewise((x, y >= x), (y, True)).replace(
        sp.Piecewise, _piecewise_to_minmax
    ) == sp.Min(x, y)
    assert sp.Piecewise((x, y < x), (y, True)).replace(
        sp.Piecewise, _piecewise_to_minmax
    ) == sp.Max(x, y)
    assert sp.Piecewise((x, y <= x), (y, True)).replace(
        sp.Piecewise, _piecewise_to_minmax
    ) == sp.Max(x, y)

    # can't replace
    assert sp.Piecewise((z, y <= x), (y, True)).replace(
        sp.Piecewise, _piecewise_to_minmax
    ) == sp.Piecewise((z, y <= x), (y, True))

    # replace recursively
    expr = sp.Piecewise(
        (sp.Piecewise((x, x < y), (y, True)), x < z),
        (sp.Piecewise((y, y < z), (z, True)), True),
    )
    replaced = expr.replace(sp.Piecewise, _piecewise_to_minmax)
    expected = sp.Piecewise(
        (sp.Min(x, y), x < z),
        (sp.Min(y, z), True),
    )
    assert replaced == expected
