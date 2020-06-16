"""
Finite Difference Check
-----------------------
This module provides functions to automatically check correctness of amici
computed sensitivities using finite difference approximations
"""

from . import (runAmiciSimulation, SensitivityOrder_none, AMICI_SUCCESS,
               SensitivityMethod_forward, Model, Solver, ExpData, ReturnData)
import numpy as np
import copy

from typing import Callable, Optional, List, Sequence


def check_finite_difference(x0: Sequence[float],
                            model: Model,
                            solver: Solver,
                            edata: ExpData,
                            ip: int,
                            fields: List[str],
                            assert_fun: Callable,
                            atol: Optional[float] = 1e-4,
                            rtol: Optional[float] = 1e-4,
                            epsilon: Optional[float] = 1e-3) -> None:
    """
    Checks the computed sensitivity based derivatives against a finite
    difference approximation.

    :param x0:
        parameter value at which to check finite difference approximation

    :param model:
        amici model

    :param solver:
        amici solver

    :param edata:
        exp data

    :param ip:
        parameter index

    :param fields:
        rdata fields for which to check the gradient

    :param assert_fun:
        function that asserts the return values of comparison, enables
        passing of custom assert function from testing frameworks

    :param atol:
        absolute tolerance for comparison

    :param rtol:
        relative tolerance for comparison

    :param epsilon:
        finite difference step-size

    """
    old_sensitivity_order = solver.getSensitivityOrder()
    old_parameters = model.getParameters()
    old_plist = model.getParameterList()

    # sensitivity
    p = copy.deepcopy(x0)
    plist = [ip]

    model.setParameters(p)
    model.setParameterList(plist)
    rdata = runAmiciSimulation(model, solver, edata)
    assert_fun(rdata['status'] == AMICI_SUCCESS)

    # finite difference
    solver.setSensitivityOrder(SensitivityOrder_none)

    # forward:
    p = copy.deepcopy(x0)
    p[ip] += epsilon/2
    model.setParameters(p)
    rdataf = runAmiciSimulation(model, solver, edata)
    assert_fun(rdata['status'] == AMICI_SUCCESS)

    # backward:
    p = copy.deepcopy(x0)
    p[ip] -= epsilon/2
    model.setParameters(p)
    rdatab = runAmiciSimulation(model, solver, edata)
    assert_fun(rdata['status'] == AMICI_SUCCESS)

    for field in fields:
        sensi_raw = rdata[f's{field}']
        fd = (rdataf[field]-rdatab[field])/epsilon
        if len(sensi_raw.shape) == 1:
            sensi = sensi_raw[0]
        elif len(sensi_raw.shape) == 2:
            sensi = sensi_raw[:, 0]
        elif len(sensi_raw.shape) == 3:
            sensi = sensi_raw[:, 0, :]
        else:
            assert_fun(False)  # not implemented
            return

        check_close(sensi, fd, assert_fun, atol, rtol, field, ip=ip)

    solver.setSensitivityOrder(old_sensitivity_order)
    model.setParameters(old_parameters)
    model.setParameterList(old_plist)


def check_derivatives(model: Model,
                      solver: Solver,
                      edata: ExpData,
                      assert_fun: Callable,
                      atol: Optional[float] = 1e-4,
                      rtol: Optional[float] = 1e-4,
                      epsilon: Optional[float] = 1e-3) -> None:
    """
    Finite differences check for likelihood gradient.

    :param model:
        amici model

    :param solver:
        amici solver

    :param edata:
        exp data

    :param assert_fun:
        function that asserts the return values of comparison, enables
        passing of custom assert function from testing frameworks

    :param atol:
        absolute tolerance for comparison

    :param rtol:
        relative tolerance for comparison

    :param epsilon:
        finite difference step-size

    """
    p = np.array(model.getParameters())

    rdata = runAmiciSimulation(model, solver, edata)
    assert_fun(rdata['status'] == AMICI_SUCCESS)

    fields = ['llh']

    leastsquares_applicable = \
        solver.getSensitivityMethod() == SensitivityMethod_forward

    if 'ssigmay' in rdata.keys() \
            and rdata['ssigmay'] is not None \
            and rdata['ssigmay'].any():
        leastsquares_applicable = False

    if leastsquares_applicable:
        fields += ['res', 'x', 'y']

        check_results(rdata, 'FIM',
                      np.dot(rdata['sres'].transpose(), rdata['sres']),
                      assert_fun,
                      1e-8, 1e-4)
        check_results(rdata, 'sllh',
                      -np.dot(rdata['res'].transpose(), rdata['sres']),
                      assert_fun,
                      1e-8, 1e-4)
    for ip in range(len(p)):
        check_finite_difference(p, model, solver, edata, ip, fields,
                                assert_fun, atol=atol, rtol=rtol,
                                epsilon=epsilon)


def check_close(result: np.array,
                expected: np.array,
                assert_fun: Callable,
                atol: float,
                rtol: float,
                field: str,
                ip: Optional[int] = None) -> None:
    """
    Compares computed values against expected values and provides rich
    output information.

    :param result:
        computed values

    :param expected:
        expected values

    :param field:
        rdata field for which the gradient is checked, only for error reporting

    :param assert_fun:
        function that asserts the return values of comparison, enables
        passing of custom assert function from testing frameworks

    :param atol:
        absolute tolerance for comparison

    :param rtol:
        relative tolerance for comparison

    :param ip:
        parameter index

    """
    close = np.isclose(result, expected, atol=atol, rtol=rtol, equal_nan=True)

    if not close.all():
        if ip is None:
            index_str = ''
            check_type = 'Regression check  '
        else:
            index_str = f'at index ip={ip} '
            check_type = 'FD check '
        print(f'{check_type} failed for {field} {index_str}for '
              f'{close.size - close.sum()} indices:')
        adev = np.abs(result - expected)
        rdev = np.abs((result - expected)/(expected + atol))
        print(f'max(adev): {adev.max()}, max(rdev): {rdev.max()}')

    assert_fun(close.all())


def check_results(rdata: ReturnData,
                  field: str,
                  expected: np.array,
                  assert_fun: Callable,
                  atol: float,
                  rtol: float) -> None:
    """
    Checks whether rdata[field] agrees with expected according to provided
    tolerances.

    :param rdata:
        simulation results as returned by
        :meth:`amici.amici.runAmiciSimulation`

    :param field:
        name of the field to check

    :param expected:
        expected values

    :param assert_fun:
        function that asserts the return values of comparison, enables
        passing of custom assert function from testing frameworks

    :param atol:
        absolute tolerance for comparison

    :param rtol:
        relative tolerance for comparison
    """

    result = rdata[field]
    if type(result) is float:
        result = np.array(result)

    check_close(result, expected, assert_fun, atol, rtol, field)
