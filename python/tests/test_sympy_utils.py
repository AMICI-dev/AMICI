"""Tests related to the sympy_utils module."""

import pytest
import sympy as sp
from amici._symbolic.sympy_utils import (
    _custom_pow_eval_derivative,
    _get_mp_context,
    _get_n_procs,
    _get_pool,
    _monkeypatched,
    _parallel_applyfunc,
    _piecewise_to_minmax,
    _shutdown_pool,
    smart_jacobian,
)
from amici.testing import skip_on_valgrind


@pytest.fixture
def nprocs(request, monkeypatch):
    """Set ``AMICI_IMPORT_NPROCS`` and clean up any worker pool afterwards."""
    monkeypatch.setenv("AMICI_IMPORT_NPROCS", str(request.param))
    yield request.param
    _shutdown_pool()


def _simplify(x):
    """Module-level (i.e. picklable) simplification function."""
    return sp.simplify(x)


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


@skip_on_valgrind
@pytest.mark.parametrize("nprocs", [1, 2], indirect=True)
def test_smart_jacobian(nprocs):
    """Serial and parallel jacobians must agree with sympy's."""
    # enough elements to exceed the serial-processing threshold
    x = sp.Matrix(sp.symbols("x0:6"))
    eq = sp.Matrix([xi**2 * xj + sp.sin(xi) for xi in x for xj in x])
    expected = eq.jacobian(x)

    actual = smart_jacobian(eq, x)

    assert actual.shape == expected.shape
    assert (sp.Matrix(actual) - expected).is_zero_matrix

    # empty input gives an empty matrix, not an error
    empty = smart_jacobian(sp.Matrix(0, 1, []), x)
    assert empty.shape == (0, len(x))


@skip_on_valgrind
@pytest.mark.parametrize("nprocs", [1, 2], indirect=True)
def test_smart_jacobian_non_symbol_vars(nprocs):
    """The sparsity pattern must also be correct for non-``Symbol`` vars."""
    t = sp.Symbol("t")
    funcs = sp.Matrix([sp.Function("x")(t), sp.Function("y")(t)])
    eq = sp.Matrix([f**2 for f in funcs] * 10)

    actual = smart_jacobian(eq, funcs)

    assert (sp.Matrix(actual) - eq.jacobian(funcs)).is_zero_matrix


@skip_on_valgrind
@pytest.mark.parametrize("nprocs", [1, 2], indirect=True)
@pytest.mark.parametrize("n_elements", [4, 40])
@pytest.mark.parametrize("sparse", [True, False])
def test_parallel_applyfunc(nprocs, n_elements, sparse):
    """``_parallel_applyfunc`` must match ``Matrix.applyfunc``.

    Tested for dense/sparse matrices and for element counts below and above
    the threshold for switching to serial processing.
    """
    x = sp.Symbol("x")
    entries = [x ** (i + 2) / x for i in range(n_elements)]
    if sparse:
        obj = sp.MutableSparseMatrix(
            n_elements,
            n_elements,
            {(i, i): entry for i, entry in enumerate(entries)},
        )
    else:
        obj = sp.MutableDenseMatrix(n_elements, 1, entries)

    actual = _parallel_applyfunc(obj, _simplify)

    assert type(actual) is type(obj)
    assert actual == obj.applyfunc(_simplify)


@skip_on_valgrind
@pytest.mark.parametrize("nprocs", [2], indirect=True)
def test_parallel_applyfunc_unpicklable(nprocs):
    """Unpicklable functions must give an actionable error message."""
    x = sp.Symbol("x")
    obj = sp.MutableDenseMatrix(40, 1, [x**2 / x] * 40)

    with pytest.raises(ValueError, match="Couldn't pickle"):
        _parallel_applyfunc(obj, lambda e: sp.simplify(e))


@skip_on_valgrind
@pytest.mark.parametrize("nprocs", [2], indirect=True)
def test_pool_is_reused(nprocs):
    """The worker pool must be created once and reused."""
    pool = _get_pool(2)
    assert _get_pool(2) is pool
    # ... but recreated if a different number of processes is requested
    assert _get_pool(3) is not pool


@skip_on_valgrind
def test_get_mp_context():
    import multiprocessing

    available = multiprocessing.get_all_start_methods()
    expected = "forkserver" if "forkserver" in available else "spawn"
    # in particular, never fork the (potentially multi-threaded) main process
    assert _get_mp_context().get_start_method() == expected


@skip_on_valgrind
def test_get_n_procs(monkeypatch):
    monkeypatch.delenv("AMICI_IMPORT_NPROCS", raising=False)
    assert _get_n_procs() == 1

    monkeypatch.setenv("AMICI_IMPORT_NPROCS", "4")
    assert _get_n_procs() == 4

    for invalid in ("0", "-1", "some_string", ""):
        monkeypatch.setenv("AMICI_IMPORT_NPROCS", invalid)
        with pytest.raises(ValueError, match="AMICI_IMPORT_NPROCS"):
            _get_n_procs()
