"""Tests for static expressions (``ws``) in :class:`DEModel`."""

from pathlib import Path

import numpy as np
import pytest
import sympy as sp
from amici._symbolic.de_model import DEModel
from amici._symbolic.de_model_components import (
    DifferentialState,
    Event,
    Expression,
    FixedParameter,
    FreeParameter,
    LogLikelihoodY,
    Observable,
    SigmaY,
    StaticExpression,
)
from amici.importers.utils import amici_time_symbol as t
from amici.sim.sundials import (
    AMICI_SUCCESS,
    ExpData,
    SensitivityMethod,
    SensitivityOrder,
    import_model_module,
    run_simulation,
)
from amici.sim.sundials.gradient_check import check_derivatives
from amici.testing import skip_on_valgrind
from numpy.testing import assert_allclose

pytestmark = pytest.mark.filterwarnings(
    # https://github.com/AMICI-dev/AMICI/issues/18
    "ignore:Adjoint sensitivity analysis for models with discontinuous "
    "right hand sides .*:UserWarning",
)

p1, p2, k1 = sp.symbols("p1 p2 k1")
q1, q2, q3 = sp.symbols("q1 q2 q3")
x1, x2 = sp.symbols("x1 x2")
e1, e2 = sp.symbols("e1 e2")
y1, my1, sigma_y1, J_y1 = sp.symbols("y1 my1 sigma_y1 J_y1")

# static expressions, deliberately listed out of dependency order:
#  q3 depends on q2 and the initial time, q2 on q1 and the fixed parameter
WS_DEFS = {
    q3: q2 * t,
    q2: 2 * q1 + k1,
    q1: p1**2,
}


def _build_model(inline_ws: bool, t0: float) -> DEModel:
    """Build a test model using static expressions in every consumer.

    With ``inline_ws=True``, the static expressions are substituted into all
    equations (with the initial time fixed to ``t0``), yielding an equivalent
    reference model without static expressions.
    """
    model = DEModel()
    model.add_component(FreeParameter(p1, "p1", 0.5))
    model.add_component(FreeParameter(p2, "p2", 0.4))
    model.add_component(FixedParameter(k1, "k1", 3.0))

    def ws(expr):
        if not inline_ws:
            return expr
        for _ in range(len(WS_DEFS)):
            expr = expr.subs(WS_DEFS)
        return expr.subs(t, t0)

    if not inline_ws:
        for sym, val in WS_DEFS.items():
            model.add_component(StaticExpression(sym, str(sym), val))

    # x0, xdot (via w and directly), y, sigmay, event trigger and
    # event assignment all reference static expressions
    # (both states have a steady state, for pre-equilibration)
    model.add_component(DifferentialState(x1, "x1", ws(q2), -e1 + ws(q3)))
    model.add_component(DifferentialState(x2, "x2", sp.Float(1.0), e1 - x2))
    model.add_component(Expression(e1, "e1", p2 * x1 * ws(q1)))
    # a static expression-dependent, but otherwise static `w` entry
    model.add_component(Expression(e2, "e2", p2 * ws(q2)))
    model.add_component(
        Observable(y1, "y1", x1 + ws(q2) + x2, measurement_symbol=my1)
    )
    model.add_component(SigmaY(sigma_y1, "sigma_y1", 0.5 + ws(q1)))
    model.add_component(
        LogLikelihoodY(
            J_y1,
            "J_y1",
            0.5 * sp.log(2 * sp.pi * sigma_y1**2)
            + 0.5 * ((y1 - my1) / sigma_y1) ** 2,
        )
    )
    model.add_component(
        Event(
            sp.Symbol("event_1"),
            "event_1",
            t - 1 - ws(q2),
            use_values_from_trigger_time=True,
            assignments={x2: x2 + 1 + ws(q3)},
        )
    )
    model.generate_basic_variables()
    return model


@skip_on_valgrind
def test_static_expression_equations():
    """Symbolic checks of the ``ws``-related equations."""
    model = _build_model(inline_ws=False, t0=1.0)

    # topologically sorted; `t` is routed through a dedicated t0 entry
    t0 = model.sym("ws")[1]
    assert t0.name == "t0"
    assert list(model.sym("ws")) == [q1, t0, q2, q3]
    assert list(model.eq("ws")) == [p1**2, t, 2 * q1 + k1, q2 * t0]

    # total derivative w.r.t. p, chained through ws
    dwsdp = model.eq("dwsdp")
    assert dwsdp.shape == (4, 2)
    assert sp.simplify(dwsdp[0, 0] - 2 * p1) == 0
    assert dwsdp[1, 0] == 0
    assert sp.simplify(dwsdp[2, 0] - 4 * p1) == 0
    assert sp.simplify(dwsdp[3, 0] - 4 * p1 * t0) == 0
    assert dwsdp[:, 1] == sp.zeros(4, 1)

    # x0 = q2 -> sx0 = dq2/dp1 = 4 p1
    assert list(model.eq("x0")) == [q2, 1.0]
    assert sp.simplify(model.eq("sx0")[0, 0] - 4 * p1) == 0

    # fixed-parameter dependence is detected through ws
    model.eq("x0_fixedParameters")
    assert model._x0_fixedParameters_idx == [0]
    assert model.state_has_fixed_parameter_initial_condition(0)
    assert not model.state_has_fixed_parameter_initial_condition(1)
    assert model.num_state_reinits() == 1

    # k-dependent ws entries
    assert list(model.eq("ws_fixedParameters")) == [2 * q1 + k1, q2 * t0]
    assert model._ws_fixedParameters_idx == [2, 3]

    # chain rule through ws in the various derivatives
    assert sp.simplify(model.eq("dwdp")[0, 0] - 2 * p1 * p2 * x1) == 0
    assert sp.simplify(model.eq("dydp")[0, 0] - 4 * p1) == 0
    assert sp.simplify(model.eq("dsigmaydp")[0, 0] - 2 * p1) == 0
    assert sp.simplify(model.eq("dxdotdp_explicit")[0, 0] - 4 * p1 * t0) == 0
    assert sp.simplify(model.eq("drootdp")[0, 0] + 4 * p1) == 0
    assert sp.simplify(model.eq("ddeltaxdp")[0][1, 0] - 4 * p1 * t0) == 0

    # `e1` depends on a state and is dynamic; `e2` only on `p` and `ws`
    assert model.static_indices("w") == [1]


def _compile(model: DEModel, name: str, tempdir) -> "amici.Model":  # noqa: F821
    from amici.exporters.sundials.de_export import DEExporter

    outdir = Path(tempdir) / name
    exporter = DEExporter(model, output_dir=outdir, model_name=name)
    exporter.generate_model_code()
    exporter.compile_model()
    return import_model_module(name, outdir).get_model()


@skip_on_valgrind
@pytest.mark.parametrize("t0", [0.0, 1.0])
def test_static_expressions_simulation(tempdir, t0):
    """Static expressions reproduce the inlined reference model, including
    sensitivities and fixed-parameter-dependent reinitialization after
    pre-equilibration."""
    model_ws = _compile(
        _build_model(inline_ws=False, t0=t0), f"ws_model_t0_{int(t0)}", tempdir
    )
    model_ref = _compile(
        _build_model(inline_ws=True, t0=t0), f"ws_ref_t0_{int(t0)}", tempdir
    )
    assert model_ws.nws == 4
    assert model_ref.nws == 0

    # avoid output timepoints coinciding with the event time (t=3), where
    #  finite differences would straddle the jump
    timepoints = np.linspace(0, 4, 9) + t0 + 0.25
    rdatas = {}
    for label, model in (("ws", model_ws), ("ref", model_ref)):
        model.set_t0(t0)
        model.set_timepoints(timepoints)
        # with this value, the event triggers inside the simulated interval
        model.set_fixed_parameters([1.5])
        model.set_reinitialize_fixed_parameter_initial_states(True)
        solver = model.create_solver()
        solver.set_relative_tolerance(1e-12)
        solver.set_absolute_tolerance(1e-14)
        # tight steady-state tolerances, so that the pre-equilibration result
        #  does not introduce noise into the finite differences
        solver.set_relative_tolerance_steady_state(1e-12)
        solver.set_absolute_tolerance_steady_state(1e-14)
        solver.set_relative_tolerance_steady_state_sensi(1e-12)
        solver.set_absolute_tolerance_steady_state_sensi(1e-14)
        solver.set_sensitivity_order(SensitivityOrder.first)
        solver.set_sensitivity_method(SensitivityMethod.forward)

        # plain simulation, to generate data
        rdata = run_simulation(model, solver)
        assert rdata.status == AMICI_SUCCESS
        # the event must actually trigger inside the simulated interval
        assert np.any(np.diff(rdata.x[:, 1]) > 0.5)
        if label == "ws":
            # `ws` is appended to `rdata.w`, with matching expression IDs
            assert list(rdata.expression_ids) == [
                "e1",
                "e2",
                "q1",
                "t0",
                "q2",
                "q3",
            ]
            assert rdata.w.shape == (len(timepoints), model.nw + model.nws)
            q1_exp = 0.5**2
            q2_exp = 2 * q1_exp + 1.5
            assert_allclose(rdata.by_id("q1"), q1_exp)
            assert_allclose(rdata.by_id("t0"), t0)
            assert_allclose(rdata.by_id("q2"), q2_exp)
            assert_allclose(rdata.by_id("q3"), q2_exp * t0)
            assert_allclose(rdata.by_id("e2"), 0.4 * q2_exp)
        edata = ExpData(rdata, 0.1, 0.0, 1)
        edata.t_start = t0
        edata.fixed_parameters = [1.5]
        # pre-equilibration with a different fixed parameter + reinit
        edata.fixed_parameters_pre_equilibration = [3.0]
        edata.reinitialize_fixed_parameter_initial_states = True

        res = {}
        for meth in (SensitivityMethod.forward, SensitivityMethod.adjoint):
            solver.set_sensitivity_method(meth)
            rdata = run_simulation(model, solver, edata)
            assert rdata.status == AMICI_SUCCESS
            res[meth] = rdata
            if label == "ws":
                # fixed-parameter-dependent entries are re-evaluated with the
                #  main simulation's fixed parameters after pre-equilibration
                assert_allclose(rdata.by_id("q1"), 0.5**2)
                assert_allclose(rdata.by_id("q2"), 2 * 0.5**2 + 1.5)
            if meth == SensitivityMethod.forward:
                # `x_ss` is only accurate to the steady-state tolerance,
                #  which makes its finite differences noisy; it is compared
                #  against the reference model below instead
                check_derivatives(
                    model,
                    solver=solver,
                    edata=edata,
                    atol=1e-5,
                    rtol=1e-4,
                    epsilon=1e-4,
                    skip_fields=["x_ss"],
                )
        rdatas[label] = res

    fsa_ws = rdatas["ws"][SensitivityMethod.forward]
    fsa_ref = rdatas["ref"][SensitivityMethod.forward]
    assert not np.allclose(fsa_ws.sllh, 0.0)
    assert_allclose(fsa_ws.x_ss, fsa_ref.x_ss, rtol=1e-6, atol=1e-8)
    assert_allclose(fsa_ws.sx_ss, fsa_ref.sx_ss, rtol=1e-6, atol=1e-8)
    assert_allclose(fsa_ws.x, fsa_ref.x, rtol=1e-8, atol=1e-10)
    assert_allclose(fsa_ws.y, fsa_ref.y, rtol=1e-8, atol=1e-10)
    assert_allclose(fsa_ws.sx, fsa_ref.sx, rtol=1e-6, atol=1e-8)
    assert_allclose(fsa_ws.sllh, fsa_ref.sllh, rtol=1e-6, atol=1e-8)
    asa_ws = rdatas["ws"][SensitivityMethod.adjoint]
    assert_allclose(asa_ws.sllh, fsa_ws.sllh, rtol=1e-5, atol=1e-7)
