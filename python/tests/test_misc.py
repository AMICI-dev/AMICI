"""Miscellaneous AMICI Python interface tests"""

import os
import subprocess

import libsbml
import pytest
import sympy as sp

import amici
from amici.ode_export import _custom_pow_eval_derivative, _monkeypatched, \
    smart_subs_dict
from amici.testing import TemporaryDirectoryWinSafe as TemporaryDirectory


def test_parameter_scaling_from_int_vector():
    """Ensure we can generate a ParameterScaling vector from Python"""

    scale_vector = amici.parameterScalingFromIntVector(
        [
            amici.ParameterScaling.log10,
            amici.ParameterScaling.ln,
            amici.ParameterScaling.none
        ])

    assert scale_vector[0] == amici.ParameterScaling.log10
    assert scale_vector[1] == amici.ParameterScaling.ln
    assert scale_vector[2] == amici.ParameterScaling.none


def test_hill_function_dwdx():
    """Kinetic laws with Hill functions, may lead to NaNs in the Jacobian
    if involved states are zero if not properly arranged symbolically.
    Test that what we are applying the right sympy simplification."""

    w = sp.Matrix([[sp.sympify('Pow(x1, p1) / (Pow(x1, p1) + a)')]])
    dwdx = w.diff(sp.Symbol('x1'))

    # Verify that without simplification we fail
    with pytest.raises(ZeroDivisionError):
        with sp.evaluate(False):
            res = dwdx.subs({'x1': 0.0})
        _ = str(res)

    # Test that powsimp does the job
    dwdx = dwdx.applyfunc(lambda x: sp.powsimp(x, deep=True))
    with sp.evaluate(False):
        res = dwdx.subs({'x1': 0.0})
    _ = str(res)


@pytest.mark.skipif(os.environ.get('AMICI_SKIP_CMAKE_TESTS', '') == 'TRUE',
                    reason='skipping cmake based test')
def test_cmake_compilation(sbml_example_presimulation_module):
    """Check that CMake build succeeds for one of the models generated during
    Python tests"""

    source_dir = os.path.dirname(sbml_example_presimulation_module.__path__[0])

    cmd = f"set -e; cd {source_dir}; mkdir -p build; cd build; "\
          "cmake ..; make"

    subprocess.run(cmd, shell=True, check=True,
                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def test_smart_subs_dict():
    expr_str = 'c + d'
    subs_dict = {
        'c': 'a + b',
        'd': 'c + a',
    }
    expected_default_str = '3*a + 2*b'
    expected_reverse_str = '2*a + b + c'

    expr_sym = sp.sympify(expr_str)
    subs_sym = {sp.sympify(k): sp.sympify(v) for k, v in subs_dict.items()}
    expected_default = sp.sympify(expected_default_str)
    expected_reverse = sp.sympify(expected_reverse_str)

    result_default = smart_subs_dict(expr_sym, subs_sym)
    result_reverse = smart_subs_dict(expr_sym, subs_sym, reverse=False)

    assert sp.simplify(result_default - expected_default).is_zero
    assert sp.simplify(result_reverse - expected_reverse).is_zero


def test_monkeypatch():
    t = sp.Symbol('t')
    n = sp.Symbol('n')
    vals = [(t, 0),
            (n, 1)]

    # check that the removable singularity still exists
    assert (t**n).diff(t).subs(vals) is sp.nan

    # check that we can monkeypatch it out
    with _monkeypatched(sp.Pow, '_eval_derivative',
                        _custom_pow_eval_derivative):
        assert (t ** n).diff(t).subs(vals) is not sp.nan

    # check that the monkeypatch is transient
    assert (t ** n).diff(t).subs(vals) is sp.nan
