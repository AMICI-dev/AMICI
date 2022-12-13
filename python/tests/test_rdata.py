"""Test amici.ReturnData(View)-related functionality"""
import numpy as np
import pytest

import amici
from numpy.testing import assert_array_equal

@pytest.fixture(scope='session')
def rdata_by_id_fixture(sbml_example_presimulation_module):
    model_module = sbml_example_presimulation_module
    model = model_module.getModel()
    model.setTimepoints(np.linspace(0, 60, 61))
    solver = model.getSolver()
    rdata = amici.runAmiciSimulation(model, solver)
    assert rdata.status == amici.AMICI_SUCCESS
    return model, rdata


def test_rdata_by_id(rdata_by_id_fixture):
    model, rdata = rdata_by_id_fixture

    assert_array_equal(
        rdata.by_id(model, 'x', model.getStateIds()[1]),
        rdata.x[:, 1]
    )

    assert_array_equal(
        rdata.by_id(model, 'y', model.getObservableIds()[0]),
        rdata.y[:, 0]
    )

    assert_array_equal(
        rdata.by_id(model, 'w', model.getExpressionIds()[1]),
        rdata.w[:, 1]
    )
