"""Adjoint sensitivity analysis for models with algebraic rules (DAEs)."""

import numpy as np
from amici import import_model_module
from amici.adapters.fiddy import run_simulation_to_function_and_derivative
from amici.importers.antimony import antimony2sbml
from amici.importers.sbml import MeasurementChannel, SbmlImporter
from amici.sim.sundials import (
    AMICI_SUCCESS,
    ExpData,
    SensitivityMethod,
    SensitivityOrder,
    run_simulation,
)
from amici.testing import skip_on_valgrind
from fiddy import check_gradient


def _check_adjoint_gradient(model):
    """Compare adjoint and finite-difference gradients of the llh."""
    model.set_timepoints(np.linspace(0, 5, 6))
    solver = model.create_solver()
    solver.set_relative_tolerance(1e-12)
    solver.set_absolute_tolerance(1e-14)

    rdata = run_simulation(model, solver)
    assert rdata.status == AMICI_SUCCESS
    edata = ExpData(rdata, 1.0, 0.0, 0)
    # perturb the data to get a non-zero gradient
    edata.set_measurements(np.asarray(edata.get_measurements()) * 1.1 + 0.05)

    solver.set_sensitivity_order(SensitivityOrder.first)
    solver.set_sensitivity_method(SensitivityMethod.adjoint)
    rdata = run_simulation(model, solver, edata)
    assert rdata.status == AMICI_SUCCESS
    # all parameters affect the objective
    assert np.all(np.abs(rdata.sllh) > 1e-3)

    function, derivative = run_simulation_to_function_and_derivative(
        amici_model=model,
        amici_solver=solver,
        amici_edata=edata,
        derivative_variables=["llh"],
    )
    point = np.asarray(model.get_free_parameters())
    result = check_gradient(function, point, derivative(point)["llh"])
    result.assert_success(always_print=True)


@skip_on_valgrind
def test_dae_adjoint_differential_observable(tempdir):
    """Adjoint sensitivities for a DAE with an observed differential state."""
    ant_model = """
    model dae_asa_diff_obs
        const k = 0.5, p = 2, q = 0.3
        x = 1; y = 2
        x' = -k * x + p * q
        0 = y - p * x
    end
    """
    module_name = "dae_asa_diff_obs"
    SbmlImporter(antimony2sbml(ant_model), from_file=False).sbml2amici(
        module_name,
        tempdir,
        observation_model=[MeasurementChannel(id_="obs_x", formula="x")],
    )
    model = import_model_module(module_name, tempdir).get_model()
    _check_adjoint_gradient(model)
