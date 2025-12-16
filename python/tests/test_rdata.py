"""Test amici.ReturnData(View)-related functionality"""

import numpy as np
import pytest
from amici.sim.sundials import (
    AMICI_SUCCESS,
    SensitivityMethod,
    SensitivityOrder,
    evaluate,
    run_simulation,
)
from amici.testing import skip_on_valgrind
from numpy.testing import assert_almost_equal, assert_array_equal


@pytest.fixture(scope="session")
def rdata_by_id_fixture(sbml_example_presimulation_module):
    model_module = sbml_example_presimulation_module
    model = model_module.get_model()
    model.set_timepoints(np.linspace(0, 60, 61))
    solver = model.create_solver()
    solver.set_sensitivity_method(SensitivityMethod.forward)
    solver.set_sensitivity_order(SensitivityOrder.first)
    rdata = run_simulation(model, solver)
    assert rdata.status == AMICI_SUCCESS
    return model, rdata


@skip_on_valgrind
def test_rdata_by_id(rdata_by_id_fixture):
    model, rdata = rdata_by_id_fixture

    assert_array_equal(rdata.by_id(model.get_state_ids()[1]), rdata.x[:, 1])
    assert_array_equal(
        rdata.by_id(model.get_state_ids()[1], "x"), rdata.x[:, 1]
    )
    assert_array_equal(
        rdata.by_id(model.get_observable_ids()[0], "y"), rdata.y[:, 0]
    )
    assert_array_equal(
        rdata.by_id(model.get_expression_ids()[1]), rdata.w[:, 1]
    )
    assert_array_equal(
        rdata.by_id(model.get_expression_ids()[1], "w"), rdata.w[:, 1]
    )
    assert_array_equal(
        rdata.by_id(model.get_state_ids()[1], "sx"), rdata.sx[:, :, 1]
    )


@skip_on_valgrind
def test_evaluate(rdata_by_id_fixture):
    # get IDs of model components
    model, rdata = rdata_by_id_fixture
    expr0_id = model.get_expression_ids()[0]
    state1_id = model.get_state_ids()[1]
    observable0_id = model.get_observable_ids()[0]

    # ensure `evaluate` works for atoms
    expr0 = rdata.by_id(expr0_id)
    assert_array_equal(expr0, evaluate(expr0_id, rdata=rdata))

    state1 = rdata.by_id(state1_id)
    assert_array_equal(state1, evaluate(state1_id, rdata=rdata))

    observable0 = rdata.by_id(observable0_id)
    assert_array_equal(observable0, evaluate(observable0_id, rdata=rdata))

    # ensure `evaluate` works for expressions
    assert_almost_equal(
        expr0 + state1 * observable0,
        evaluate(f"{expr0_id} + {state1_id} * {observable0_id}", rdata=rdata),
    )
