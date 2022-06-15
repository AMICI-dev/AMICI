"""Tests related to amici.ExpData via Python"""
import numpy as np

import amici

from test_sbml_import import model_units_module


def test_edata_sensi_unscaling(model_units_module):
    """
    ExpData parameters should be used for unscaling initial state
    sensitivities.
    """
    parameters0 = (5, 5)
    parameters1 = (2, 2)

    sx0 = np.array((3, 3, 3, 3))

    parameter_scales_log10 = \
        [amici.ParameterScaling.log10.value]*len(parameters0)
    amici_parameter_scales_log10 = \
        amici.parameterScalingFromIntVector(parameter_scales_log10)

    model = model_units_module.getModel()
    model.setTimepoints(np.linspace(0, 1, 3))
    model.setParameterScale(parameter_scales_log10)
    model.setParameters(parameters0)

    solver = model.getSolver()
    solver.setSensitivityOrder(amici.SensitivityOrder.first)

    edata = amici.ExpData(model)
    edata.pscale = amici_parameter_scales_log10
    edata.parameters = parameters0
    edata.sx0 = sx0

    edata2 = amici.ExpData(model)
    edata2.pscale = amici_parameter_scales_log10
    edata2.parameters = parameters1
    edata2.sx0 = sx0

    rdata = amici.runAmiciSimulation(model, solver, edata)
    rdata2 = amici.runAmiciSimulation(model, solver, edata2)

    # The initial state sensitivities are as specified.
    assert (rdata.sx0.flatten() == sx0).all()
    assert (rdata2.sx0.flatten() == sx0).all()
