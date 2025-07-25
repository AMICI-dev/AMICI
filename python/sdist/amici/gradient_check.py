"""
Finite Difference Check
-----------------------
This module provides functions to automatically check correctness of amici
computed sensitivities using finite difference approximations
"""

import copy
from collections.abc import Sequence

import numpy as np

from . import (
    AMICI_SUCCESS,
    ExpData,
    Model,
    ParameterScaling,
    ReturnData,
    SensitivityMethod,
    SensitivityOrder,
    Solver,
    runAmiciSimulation,
)


def check_finite_difference(
    x0: Sequence[float],
    model: Model,
    solver: Solver,
    edata: ExpData,
    ip: int,
    fields: list[str],
    atol: float | None = 1e-4,
    rtol: float | None = 1e-4,
    epsilon: float | None = 1e-3,
) -> None:
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

    :param atol:
        absolute tolerance for comparison

    :param rtol:
        relative tolerance for comparison

    :param epsilon:
        finite difference step-size

    """
    og_sensitivity_order = solver.getSensitivityOrder()
    og_parameters = model.getParameters()
    og_plist = model.getParameterList()
    if edata:
        og_eplist = edata.plist

    # sensitivity
    p = copy.deepcopy(x0)
    plist = [ip]

    model.setParameters(p)
    model.setParameterList(plist)
    if edata:
        edata.plist = plist

    # simulation with gradient
    if int(og_sensitivity_order) < int(SensitivityOrder.first):
        solver.setSensitivityOrder(SensitivityOrder.first)
    rdata = runAmiciSimulation(model, solver, edata)
    if rdata["status"] != AMICI_SUCCESS:
        raise AssertionError(f"Simulation failed (status {rdata['status']}")

    # finite difference
    solver.setSensitivityOrder(SensitivityOrder.none)

    pf = copy.deepcopy(x0)
    pb = copy.deepcopy(x0)
    pscale = model.getParameterScale()[ip]
    if x0[ip] == 0 or pscale != int(ParameterScaling.none):
        pf[ip] += epsilon / 2
        pb[ip] -= epsilon / 2
    else:
        pf[ip] *= 1 + epsilon / 2
        pb[ip] /= 1 + epsilon / 2

    # forward:
    model.setParameters(pf)
    rdataf = runAmiciSimulation(model, solver, edata)
    if rdataf["status"] != AMICI_SUCCESS:
        raise AssertionError(f"Simulation failed (status {rdataf['status']}")

    # backward:
    model.setParameters(pb)
    rdatab = runAmiciSimulation(model, solver, edata)
    if rdatab["status"] != AMICI_SUCCESS:
        raise AssertionError(f"Simulation failed (status {rdatab['status']}")

    for field in fields:
        sensi_raw = rdata[f"s{field}"]
        fd = (rdataf[field] - rdatab[field]) / (pf[ip] - pb[ip])
        if len(sensi_raw.shape) == 1 or field == "x_ss":
            sensi = sensi_raw[0]
        elif len(sensi_raw.shape) == 2:
            sensi = sensi_raw[:, 0]
        elif len(sensi_raw.shape) == 3:
            sensi = sensi_raw[:, 0, :]
        else:
            raise NotImplementedError()

        _check_close(
            sensi,
            fd,
            atol=atol,
            rtol=rtol,
            field=field,
            ip=ip,
            parameter_id=model.getParameterIds()[ip]
            if model.hasParameterIds()
            else None,
        )

    solver.setSensitivityOrder(og_sensitivity_order)
    model.setParameters(og_parameters)
    model.setParameterList(og_plist)
    if edata:
        edata.plist = og_eplist


def check_derivatives(
    model: Model,
    solver: Solver,
    edata: ExpData | None = None,
    atol: float | None = 1e-4,
    rtol: float | None = 1e-4,
    epsilon: float | None = 1e-3,
    check_least_squares: bool = True,
    skip_zero_pars: bool = False,
    skip_fields: list[str] | None = None,
) -> None:
    """
    Finite differences check for likelihood gradient.

    :param model: amici model
    :param solver: amici solver
    :param edata: exp data
    :param atol: absolute tolerance for comparison
    :param rtol: relative tolerance for comparison
    :param epsilon: finite difference step-size
    :param check_least_squares: whether to check least squares related values.
    :param skip_zero_pars: whether to perform FD checks for parameters that
        are zero
    :param skip_fields: list of fields to skip
    """
    p = np.array(model.getParameters())

    og_sens_order = solver.getSensitivityOrder()

    if int(og_sens_order) < int(SensitivityOrder.first):
        solver.setSensitivityOrder(SensitivityOrder.first)
    rdata = runAmiciSimulation(model, solver, edata)
    solver.setSensitivityOrder(og_sens_order)

    if rdata["status"] != AMICI_SUCCESS:
        raise AssertionError(f"Simulation failed (status {rdata['status']}")

    fields = []

    if (
        solver.getSensitivityMethod() == SensitivityMethod.forward
        and solver.getSensitivityOrder() <= SensitivityOrder.first
    ):
        if rdata.sx_ss is not None:
            fields.append("x_ss")

        fields.append("x")

    leastsquares_applicable = (
        solver.getSensitivityMethod() == SensitivityMethod.forward
        and edata is not None
    )

    if (
        "ssigmay" in rdata.keys()
        and rdata["ssigmay"] is not None
        and rdata["ssigmay"].any()
        and not model.getAddSigmaResiduals()
    ):
        leastsquares_applicable = False

    if check_least_squares and leastsquares_applicable:
        fields += ["res", "y"]

        _check_results(
            rdata,
            "FIM",
            np.dot(rdata["sres"].T, rdata["sres"]),
            atol=1e-8,
            rtol=1e-4,
        )
        _check_results(
            rdata,
            "sllh",
            -np.dot(rdata["res"].T, rdata["sres"]),
            atol=1e-8,
            rtol=1e-4,
        )

    if edata is not None:
        fields.append("llh")

    fields = [f for f in fields if f not in (skip_fields or [])]

    # only check the sensitivities w.r.t. the selected parameters
    plist = model.getParameterList()
    if edata and edata.plist:
        plist = edata.plist

    for ip, pval in enumerate(p):
        if plist and ip not in plist:
            continue
        if pval == 0.0 and skip_zero_pars:
            continue
        check_finite_difference(
            p,
            model,
            solver,
            edata,
            ip,
            fields,
            atol=atol,
            rtol=rtol,
            epsilon=epsilon,
        )


def _check_close(
    result: np.array,
    expected: np.array,
    atol: float,
    rtol: float,
    field: str,
    ip: int | None = None,
    parameter_id: str | None = None,
    verbose: bool | None = True,
) -> None:
    """
    Compares computed values against expected values and provides rich
    output information.

    :param result:
        computed values

    :param expected:
        expected values

    :param field:
        rdata field for which the gradient is checked, only for error reporting

    :param atol:
        absolute tolerance for comparison

    :param rtol:
        relative tolerance for comparison

    :param ip:
        parameter index, for more informative output

    :param parameter_id:
        parameter ID, for more informative output

    :param verbose:
        produce a more verbose error message in case of unmatched expectations
    """
    close = np.isclose(result, expected, atol=atol, rtol=rtol, equal_nan=True)
    if close.all():
        return

    if ip is None:
        index_str = ""
        check_type = "Regression check"
    else:
        index_str = f"at index ip={ip} "
        if parameter_id:
            index_str += f"({parameter_id}) "
        check_type = "FD check"

    lines = [
        f"{check_type} failed for {field} {index_str}for "
        f"{close.size - close.sum()} indices:"
    ]
    if verbose:
        for idx in np.argwhere(~close):
            idx = tuple(idx)
            if result.shape:
                rr = result[idx]
            else:
                rr = result
            lines.append(f"\tat {idx}: Expected {expected[idx]}, got {rr}")
    adev = np.abs(result - expected)
    rdev = np.abs((result - expected) / (expected + atol))
    lines.append(f"max(adev): {adev.max()}, max(rdev): {rdev.max()}")

    raise AssertionError("\n".join(lines))


def _check_results(
    rdata: ReturnData, field: str, expected: np.array, atol: float, rtol: float
) -> None:
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

    :param atol:
        absolute tolerance for comparison

    :param rtol:
        relative tolerance for comparison
    """
    if field in ("rdrm", "messages"):
        return

    result = rdata[field]

    if isinstance(result, str):
        if result != expected:
            raise AssertionError(
                f"Expected {expected} but got {result} for field {field}"
            )
        return

    if type(result) is float:  # noqa E721
        result = np.array(result)

    _check_close(
        result=result, expected=expected, atol=atol, rtol=rtol, field=field
    )
