"""Functionality for working with sympy objects."""

import atexit
import contextlib
import logging
import os
import threading
from collections.abc import Callable
from functools import wraps
from itertools import starmap
from typing import TYPE_CHECKING, Any

import sympy as sp

from amici.logging import get_logger, log_execution_time

if TYPE_CHECKING:
    from multiprocessing.pool import Pool

logger = get_logger(__name__, logging.ERROR)

__all__ = [
    "smart_jacobian",
    "smart_multiply",
    "smart_is_zero_matrix",
    "_monkeypatch_sympy",
    "_parallel_applyfunc",
    "_piecewise_to_minmax",
]

# Number of matrix elements below which the inter-process communication
#  overhead is expected to outweigh any speed-up from parallel processing.
_MIN_PARALLEL_ELEMENTS = 16

# The worker pool for parallel model import, the ``(n_procs, pid)`` it was
#  created for, and the lock guarding both. Creating a worker pool is expensive
#  (every worker has to import sympy and amici), and model import performs many
#  parallelizable operations. Therefore, a single pool is created lazily and
#  reused for all of them, instead of creating a new one for each operation.
_pool: "Pool | None" = None
_pool_config: tuple[int, int] | None = None
_pool_lock = threading.RLock()


def _get_n_procs() -> int:
    """Number of processes to be used for model import.

    Controlled via the ``AMICI_IMPORT_NPROCS`` environment variable.
    """
    val = os.environ.get("AMICI_IMPORT_NPROCS", "1")
    try:
        n_procs = int(val)
        if n_procs < 1:
            raise ValueError
    except ValueError:
        raise ValueError(
            f"Invalid value for AMICI_IMPORT_NPROCS: {val!r}. "
            "Must be a positive integer."
        ) from None
    return n_procs


def _get_mp_context():
    """Get the multiprocessing context to be used for model import.

    ``fork`` is not used, since forking a potentially multi-threaded process
    may deadlock (see e.g. https://stackoverflow.com/a/66113051).

    ``forkserver`` is used where available, because its workers are forked from
    a small, single-threaded server process that imports sympy and amici only
    once, whereas every ``spawn`` worker has to import them from scratch. This
    makes creating a pool several times cheaper, without being affected by the
    problems of forking the main process. Where ``forkserver`` is unavailable
    (Windows), ``spawn`` is used.
    """
    from multiprocessing import get_all_start_methods, get_context

    if "forkserver" not in get_all_start_methods():
        return get_context("spawn")

    ctx = get_context("forkserver")
    # import the modules required by the workers once in the forkserver
    #  process, instead of once in every worker process
    #  (unimportable modules are silently ignored by the forkserver)
    ctx.set_forkserver_preload(["sympy", __name__])
    return ctx


def _get_pool(n_procs: int) -> "Pool":
    """Get the (lazily created) worker pool with ``n_procs`` processes."""
    global _pool, _pool_config

    with _pool_lock:
        if _pool is not None and _pool_config != (n_procs, os.getpid()):
            # either the requested number of processes changed, or we are in a
            #  forked child process that inherited the parent's unusable pool
            _shutdown_pool()

        if _pool is None:
            _pool = _get_mp_context().Pool(n_procs)
            _pool_config = (n_procs, os.getpid())

        return _pool


def _shutdown_pool() -> None:
    """Shut down the worker pool, if any."""
    global _pool, _pool_config

    with _pool_lock:
        if _pool is None:
            return
        if _pool_config[1] == os.getpid():
            # only the process that created the pool may terminate it
            _pool.terminate()
            _pool.join()
        _pool = _pool_config = None


atexit.register(_shutdown_pool)


def _custom_pow_eval_derivative(self, s):
    """
    Custom Pow derivative that removes a removable singularity for
    ``self.base == 0`` and ``self.base.diff(s) == 0``. This function is
    intended to be monkeypatched into :py:method:`sympy.Pow._eval_derivative`.

    :param self:
        sp.Pow class

    :param s:
        variable with respect to which the derivative will be computed
    """
    dbase = self.base.diff(s)
    dexp = self.exp.diff(s)
    part1 = sp.Pow(self.base, self.exp - 1) * self.exp * dbase
    part2 = self * dexp * sp.log(self.base)
    if self.base.is_nonzero or dbase.is_nonzero or part2.is_zero:
        # first piece never applies or is zero anyway
        return part1 + part2

    return part1 + sp.Piecewise(
        (self.base, sp.And(sp.Eq(self.base, 0), sp.Eq(dbase, 0))),
        (part2, True),
    )


@contextlib.contextmanager
def _monkeypatched(obj: object, name: str, patch: Any):
    """
    Temporarily monkeypatches an object.

    :param obj:
        object to be patched

    :param name:
        name of the attribute to be patched

    :param patch:
        patched value
    """
    pre_patched_value = getattr(obj, name)
    setattr(obj, name, patch)
    try:
        yield object
    finally:
        setattr(obj, name, pre_patched_value)


def _monkeypatch_sympy(func):
    """
    Decorator that temporarily monkeypatches sympy.Pow._eval_derivative.
    """

    @wraps(func)
    def wrapper(*args, **kwargs):
        with _monkeypatched(
            sp.Pow, "_eval_derivative", _custom_pow_eval_derivative
        ):
            return func(*args, **kwargs)

    return wrapper


@log_execution_time("running smart_jacobian", logger)
def smart_jacobian(
    eq: sp.MutableDenseMatrix, sym_var: sp.MutableDenseMatrix
) -> sp.MutableSparseMatrix:
    """
    Wrapper around symbolic jacobian with some additional checks that reduce
    computation time for large matrices

    :param eq:
        equation
    :param sym_var:
        differentiation variable
    :return:
        jacobian of eq wrt sym_var
    """
    nrow = eq.shape[0]
    ncol = sym_var.shape[0]
    if (
        not min(eq.shape)
        or not min(sym_var.shape)
        or smart_is_zero_matrix(eq)
        or smart_is_zero_matrix(sym_var)
    ):
        return sp.MutableSparseMatrix(nrow, ncol, dict())

    # preprocess sparsity pattern
    if all(b.is_Symbol for b in sym_var):
        # `Basic.has` traverses the full expression tree on every call.
        #  Collecting the symbols of each row once is equivalent, but avoids
        #  re-traversing every row for each of the `ncol` variables.
        symbols_by_row = [a.atoms(sp.Symbol) for a in eq]
        elements = (
            (i, j, a, b)
            for i, a in enumerate(eq)
            for j, b in enumerate(sym_var)
            if b in symbols_by_row[i]
        )
    else:
        elements = (
            (i, j, a, b)
            for i, a in enumerate(eq)
            for j, b in enumerate(sym_var)
            if a.has(b)
        )

    if (n_procs := _get_n_procs()) == 1:
        # serial
        return sp.MutableSparseMatrix(
            nrow, ncol, dict(starmap(_jacobian_element, elements))
        )

    # parallel -- the pool consumes the full iterable anyway
    elements = list(elements)
    if len(elements) < _MIN_PARALLEL_ELEMENTS:
        return sp.MutableSparseMatrix(
            nrow, ncol, dict(starmap(_jacobian_element, elements))
        )

    mapped = _get_pool(n_procs).starmap(_jacobian_element, elements)
    return sp.MutableSparseMatrix(nrow, ncol, dict(mapped))


@log_execution_time("running smart_multiply", logger)
def smart_multiply(
    x: sp.MutableDenseMatrix | sp.MutableSparseMatrix,
    y: sp.MutableDenseMatrix,
) -> sp.MutableDenseMatrix | sp.MutableSparseMatrix:
    """
    Wrapper around symbolic multiplication with some additional checks that
    reduce computation time for large matrices

    :param x:
        educt 1
    :param y:
        educt 2
    :return:
        product
    """
    if (
        not x.shape[0]
        or not y.shape[1]
        or smart_is_zero_matrix(x)
        or smart_is_zero_matrix(y)
    ):
        return sp.zeros(x.shape[0], y.shape[1])
    return x.multiply(y)


def smart_is_zero_matrix(
    x: sp.MutableDenseMatrix | sp.MutableSparseMatrix,
) -> bool:
    """A faster implementation of sympy's is_zero_matrix

    Avoids repeated indexer type checks and double iteration to distinguish
    False/None. Found to be about 100x faster for large matrices.

    :param x: Matrix to check
    """

    if isinstance(x, sp.MutableDenseMatrix):
        return all(xx.is_zero is True for xx in x.flat())

    if isinstance(x, list):
        return all(smart_is_zero_matrix(xx) for xx in x)

    return x.nnz() == 0


def _jacobian_element(i, j, eq_i, sym_var_j):
    """Compute a single element of a jacobian"""
    return (i, j), eq_i.diff(sym_var_j)


def _parallel_applyfunc(obj: sp.Matrix, func: Callable) -> sp.Matrix:
    """Parallel implementation of sympy's Matrix.applyfunc"""
    if (n_procs := _get_n_procs()) == 1:
        # serial
        return obj.applyfunc(func)

    from multiprocessing.reduction import ForkingPickler

    from sympy.matrices.dense import DenseMatrix

    if isinstance(obj, DenseMatrix):
        values = list(obj)
        dok = None
    elif isinstance(obj, sp.SparseMatrix):
        dok = obj.todok()
        values = list(dok.values())
    else:
        raise ValueError(f"Unsupported matrix type {type(obj)}")

    if len(values) < _MIN_PARALLEL_ELEMENTS:
        # not worth the inter-process communication overhead
        return obj.applyfunc(func)

    try:
        # check upfront -- passing `func` to the pool would surface this as an
        #  opaque error from the pool's task handler thread, with an exception
        #  type that depends on why exactly `func` is unpicklable
        ForkingPickler.dumps(func)
    except Exception as e:
        raise ValueError(
            f"Couldn't pickle {func}. This is likely because the argument "
            "was not a module-level function. Either rewrite the argument "
            "to a module-level function or disable parallelization by "
            "setting `AMICI_IMPORT_NPROCS=1`."
        ) from e

    mapped = _get_pool(n_procs).map(func, values)

    if dok is None:
        return obj._new(obj.rows, obj.cols, mapped)

    dok = {k: v for k, v in zip(dok.keys(), mapped, strict=True) if v != 0}
    return obj._new(obj.rows, obj.cols, dok)


def _piecewise_to_minmax(
    *expr_cond_pairs: tuple[tuple[sp.Basic, sp.Basic], ...],
) -> sp.Basic:
    """Replace min/max defined via Piecewise with plain Min/Max.

    To be used in ``expr = expr.replace(sp.Piecewise, pw_to_minmax)``.
    """
    if len(expr_cond_pairs) == 2 and expr_cond_pairs[-1][1] == sp.true:
        (expr1, cond1), (expr2, cond2) = expr_cond_pairs
        if cond1.args == (expr1, expr2) and cond1.func in (sp.Lt, sp.Le):
            return sp.Min(expr1, expr2)
        elif cond1.args == (expr1, expr2) and cond1.func in (sp.Gt, sp.Ge):
            return sp.Max(expr1, expr2)
    return sp.Piecewise(*expr_cond_pairs)
