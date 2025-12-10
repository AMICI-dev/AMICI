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
    ReturnDataView,
    SensitivityMethod,
    SensitivityOrder,
    Solver,
    run_simulation,
)

__all__ = ["check_finite_difference"]


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
    p = copy.deepcopy(x0)
    plist = [ip]

    # store original settings and apply new ones
    og_sensitivity_order = solver.get_sensitivity_order()
    og_parameters = model.get_free_parameters()
    og_plist = model.get_parameter_list()
    if edata:
        og_eplist = edata.plist
        og_eparameters = edata.free_parameters

        edata.plist = plist
        # we always set parameters via the model below
        edata.free_parameters = []
        pscale = (
            edata.pscale if len(edata.pscale) else model.get_parameter_scale()
        )
    else:
        pscale = model.get_parameter_scale()
        model.set_parameter_list(plist)

    model.set_parameter_scale(pscale)
    model.set_free_parameters(p)

    # simulation with gradient
    if int(og_sensitivity_order) < int(SensitivityOrder.first):
        solver.set_sensitivity_order(SensitivityOrder.first)
    rdata = run_simulation(model, solver, edata)
    if rdata["status"] != AMICI_SUCCESS:
        raise AssertionError(f"Simulation failed (status {rdata['status']}")

    # finite difference
    solver.set_sensitivity_order(SensitivityOrder.none)

    pf = copy.deepcopy(x0)
    pb = copy.deepcopy(x0)
    if x0[ip] == 0 or pscale[ip] != int(ParameterScaling.none):
        pf[ip] += epsilon / 2
        pb[ip] -= epsilon / 2
    else:
        pf[ip] *= 1 + epsilon / 2
        pb[ip] /= 1 + epsilon / 2

    # forward:
    model.set_free_parameters(pf)
    rdataf = run_simulation(model, solver, edata)
    if rdataf["status"] != AMICI_SUCCESS:
        raise AssertionError(f"Simulation failed (status {rdataf['status']}")

    # backward:
    model.set_free_parameters(pb)
    rdatab = run_simulation(model, solver, edata)
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

        try:
            _check_close(
                sensi,
                fd,
                atol=atol,
                rtol=rtol,
                field=field,
                ip=ip,
                parameter_id=model.get_free_parameter_ids()[ip]
                if model.has_free_parameter_ids()
                else None,
            )
        except Exception as e:
            sm = SensitivityMethod(solver.get_sensitivity_method())
            e.add_note(f"Sensitivity method was {sm!r}")
            raise e

    solver.set_sensitivity_order(og_sensitivity_order)
    model.set_free_parameters(og_parameters)
    model.set_parameter_list(og_plist)
    if edata:
        edata.plist = og_eplist
        edata.free_parameters = og_eparameters


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
    :param edata: ExpData instance. If provided, ExpData settings will
        override model settings where applicable (`plist`, `parmeters`, ...).
    :param atol: absolute tolerance for comparison
    :param rtol: relative tolerance for comparison
    :param epsilon: finite difference step-size
    :param check_least_squares: whether to check least squares related values.
    :param skip_zero_pars: whether to perform FD checks for parameters that
        are zero
    :param skip_fields: list of fields to skip
    """
    if edata and edata.free_parameters:
        p = np.array(edata.free_parameters)
    else:
        p = np.array(model.get_free_parameters())

    og_sens_order = solver.get_sensitivity_order()

    if int(og_sens_order) < int(SensitivityOrder.first):
        solver.set_sensitivity_order(SensitivityOrder.first)
    rdata = run_simulation(model, solver, edata)
    solver.set_sensitivity_order(og_sens_order)

    if rdata["status"] != AMICI_SUCCESS:
        raise AssertionError(f"Simulation failed (status {rdata['status']}")

    fields = []

    if (
        solver.get_sensitivity_method() == SensitivityMethod.forward
        and solver.get_sensitivity_order() <= SensitivityOrder.first
    ):
        if rdata.sx_ss is not None:
            fields.append("x_ss")

        fields.append("x")

    leastsquares_applicable = (
        solver.get_sensitivity_method() == SensitivityMethod.forward
        and edata is not None
    )

    if (
        "ssigmay" in rdata.keys()
        and rdata["ssigmay"] is not None
        and rdata["ssigmay"].any()
        and not model.get_add_sigma_residuals()
    ):
        leastsquares_applicable = False

    if check_least_squares and leastsquares_applicable:
        fields += ["y", "res"]

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
    plist = model.get_parameter_list()
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
    result: np.ndarray,
    expected: np.ndarray,
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
    rdata: ReturnDataView,
    field: str,
    expected: np.ndarray,
    atol: float,
    rtol: float,
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
