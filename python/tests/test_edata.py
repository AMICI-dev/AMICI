"""Tests related to amici.ExpData via Python"""

import numpy as np
import pytest
from amici.sim.sundials import (
    ExpData,
    ParameterScaling,
    SensitivityOrder,
    parameter_scaling_from_int_vector,
    run_simulation,
)
from amici.testing import skip_on_valgrind
from numpy.testing import assert_allclose


@skip_on_valgrind
@pytest.mark.usefixtures("model_units_module")
def test_edata_sensi_unscaling(model_units_module):  # noqa: F811
    """
    ExpData parameters should be used for unscaling initial state
    sensitivities.
    """
    parameters0 = (5, 5)
    parameters1 = (2, 2)

    sx0 = (3, 3, 3, 3)

    parameter_scales_log10 = [ParameterScaling.log10.value] * len(parameters0)
    amici_parameter_scales_log10 = parameter_scaling_from_int_vector(
        parameter_scales_log10
    )

    model = model_units_module.get_model()
    model.set_timepoints(np.linspace(0, 1, 3))
    model.set_parameter_scale(parameter_scales_log10)
    model.set_free_parameters(parameters0)

    solver = model.create_solver()
    solver.set_sensitivity_order(SensitivityOrder.first)

    edata0 = ExpData(model)
    edata0.pscale = amici_parameter_scales_log10
    edata0.free_parameters = parameters0
    edata0.sx0 = sx0

    edata1 = ExpData(model)
    edata1.pscale = amici_parameter_scales_log10
    edata1.free_parameters = parameters1
    edata1.sx0 = sx0

    rdata0 = run_simulation(model, solver, edata0)
    rdata1 = run_simulation(model, solver, edata1)

    # The initial state sensitivities are as specified.
    assert_allclose(rdata0.sx0.flatten(), sx0)
    assert_allclose(rdata1.sx0.flatten(), sx0)
