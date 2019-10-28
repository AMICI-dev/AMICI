from . import (runAmiciSimulation, SensitivityOrder_none,
               SensitivityMethod_forward)
import numpy as np
import copy


def check_finite_difference(x0, model, solver, edata, ip, fields,
                            assert_fun, atol=1e-4, rtol=1e-4, epsilon=1e-3):
    old_sensitivity_order = solver.getSensitivityOrder()
    old_parameters = model.getParameters()
    old_plist = model.getParameterList()

    # sensitivity
    p = copy.deepcopy(x0)
    plist = [ip]

    model.setParameters(p)
    model.setParameterList(plist)
    rdata = runAmiciSimulation(model, solver, edata)

    # finite difference
    solver.setSensitivityOrder(SensitivityOrder_none)

    # forward:
    p = copy.deepcopy(x0)
    p[ip] += epsilon/2
    model.setParameters(p)
    rdataf = runAmiciSimulation(model, solver, edata)

    # backward:
    p = copy.deepcopy(x0)
    p[ip] -= epsilon/2
    model.setParameters(p)
    rdatab = runAmiciSimulation(model, solver, edata)

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

        check_close(sensi, fd, assert_fun, atol, rtol, field, ip=ip)

    solver.setSensitivityOrder(old_sensitivity_order)
    model.setParameters(old_parameters)
    model.setParameterList(old_plist)


def check_derivatives(model, solver, edata, assert_fun,
                      atol=1e-4, rtol=1e-4, epsilon=1e-3):
    """Finite differences check for likelihood gradient

    Arguments:
        model: amici model
        solver: amici solver
        edata: exp data
        atol: absolute tolerance
        rtol: relative tolerance
        epsilon: finite difference step-size
    """
    from scipy.optimize import check_grad

    p = np.array(model.getParameters())

    rdata = runAmiciSimulation(model, solver, edata)

    fields = ['llh']

    leastsquares_applicable = \
        solver.getSensitivityMethod() == SensitivityMethod_forward

    if 'ssigmay' in rdata.keys():
        if rdata['ssigmay'] is not None:
            if rdata['ssigmay'].any():
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


def check_close(result, expected, assert_fun, atol, rtol, field, ip=None):
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


def check_results(rdata, field, expected, assert_fun, atol, rtol):
    """
    checks whether rdata[field] agrees with expected according to provided
    tolerances

    Arguments:
        rdata: simulation results as returned by amici.runAmiciSimulation
        field: name of the field to check
        expected: expected test results
        atol: absolute tolerance
        rtol: relative tolerance
    """

    result = rdata[field]
    if type(result) is float:
        result = np.array(result)

    check_close(result, expected, assert_fun, atol, rtol, field)

