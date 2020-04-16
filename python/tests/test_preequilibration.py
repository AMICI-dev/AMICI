"""Tests for preequilibration"""

import itertools

import amici
import numpy as np
import pytest
from test_pysb import get_data


@pytest.fixture
def preeq_fixture(pysb_example_presimulation_module):
    model = pysb_example_presimulation_module.getModel()
    model.setReinitializeFixedParameterInitialStates(True)

    solver = model.getSolver()
    solver.setSensitivityOrder(amici.SensitivityOrder_first)
    solver.setSensitivityMethod(amici.SensitivityMethod_forward)

    edata = get_data(model)
    edata.t_presim = 2
    edata.fixedParameters = [10, 2]
    edata.fixedParametersPresimulation = [3, 2]
    edata.fixedParametersPreequilibration = [3, 0]
    edata.setTimepoints([1, 5])

    edata_preeq = amici.ExpData(edata)
    edata_preeq.t_presim = 0
    edata_preeq.setTimepoints([np.infty])
    edata_preeq.fixedParameters = \
        edata.fixedParametersPreequilibration
    edata_preeq.fixedParametersPresimulation = ()
    edata_preeq.fixedParametersPreequilibration = ()

    edata_presim = amici.ExpData(edata)
    edata_presim.t_presim = 0
    edata_presim.setTimepoints([edata.t_presim])
    edata_presim.fixedParameters = \
        edata.fixedParametersPresimulation
    edata_presim.fixedParametersPresimulation = ()
    edata_presim.fixedParametersPreequilibration = ()

    edata_sim = amici.ExpData(edata)
    edata_sim.t_presim = 0
    edata_sim.setTimepoints(edata.getTimepoints())
    edata_sim.fixedParameters = \
        edata.fixedParameters
    edata_sim.fixedParametersPresimulation = ()
    edata_sim.fixedParametersPreequilibration = ()

    pscales = [
        amici.ParameterScaling_log10, amici.ParameterScaling_ln,
        amici.ParameterScaling_none,
        amici.parameterScalingFromIntVector([
            amici.ParameterScaling_log10, amici.ParameterScaling_ln,
            amici.ParameterScaling_none, amici.ParameterScaling_log10,
            amici.ParameterScaling_ln, amici.ParameterScaling_none
        ])
    ]

    plists = [
        [3, 1, 2, 4], [0, 1, 2, 3, 4, 5], [5, 3, 2, 0, 4, 1],
        [1, 2, 3, 4, 5], [1, 1, 1],
    ]

    return (model, solver, edata, edata_preeq,
            edata_presim, edata_sim, pscales, plists)


def test_manual_preequilibration(preeq_fixture):
    """Manual preequilibration"""

    model, solver, edata, edata_preeq, \
        edata_presim, edata_sim, pscales, plists = preeq_fixture

    settings = itertools.product(pscales, plists)

    for pscale, plist in settings:

        model.setInitialStates([])
        model.setInitialStateSensitivities([])
        model.setParameterList(plist)
        model.setParameterScale(pscale)

        # combined
        rdata_auto = amici.runAmiciSimulation(model, solver, edata)
        assert rdata_auto.status == amici.AMICI_SUCCESS

        # manual preequilibration
        rdata_preeq = amici.runAmiciSimulation(model, solver, edata_preeq)
        assert rdata_preeq.status == amici.AMICI_SUCCESS

        # manual reinitialization + presimulation
        x0 = rdata_preeq['x'][0, :]
        x0[1] = edata_presim.fixedParameters[0]
        x0[2] = edata_presim.fixedParameters[1]
        sx0 = rdata_preeq['sx'][0, :, :]
        sx0[:, 1] = 0
        sx0[:, 2] = 0
        model.setInitialStates(x0)
        model.setInitialStateSensitivities(sx0.flatten())
        rdata_presim = amici.runAmiciSimulation(model, solver, edata_presim)
        assert rdata_presim.status == amici.AMICI_SUCCESS

        # manual reinitialization + simulation
        x0 = rdata_presim['x'][0, :]
        x0[1] = edata_sim.fixedParameters[0]
        x0[2] = edata_sim.fixedParameters[1]
        sx0 = rdata_presim['sx'][0, :, :]
        sx0[:, 1] = 0
        sx0[:, 2] = 0
        model.setInitialStates(x0)
        model.setInitialStateSensitivities(sx0.flatten())
        rdata_sim = amici.runAmiciSimulation(model, solver, edata_sim)
        assert rdata_sim.status == amici.AMICI_SUCCESS

        for variable in ['x', 'sx']:
            assert np.isclose(
                rdata_auto[variable],
                rdata_sim[variable],
                1e-6, 1e-6
            ).all(), dict(pscale=pscale, plist=plist, variable=variable)


def test_parameter_reordering(preeq_fixture):
    """Test parameter reordering"""

    model, solver, edata, edata_preeq, \
        edata_presim, edata_sim, pscales, plists = preeq_fixture

    rdata_ordered = amici.runAmiciSimulation(model, solver, edata)

    for plist in plists:
        model.setParameterList(plist)
        rdata_reordered = amici.runAmiciSimulation(model, solver, edata)

        for ip, p_index in enumerate(plist):
            assert np.isclose(
                rdata_ordered['sx'][:, p_index, :],
                rdata_reordered['sx'][:, ip, :],
                1e-6, 1e-6
            ).all(), plist


def test_parameter_in_expdata(preeq_fixture):
    """Test parameter in ExpData"""

    model, solver, edata, edata_preeq, edata_presim, \
        edata_sim, pscales, plists = preeq_fixture

    rdata = amici.runAmiciSimulation(model, solver, edata)

    # get initial states will compute initial states if nothing is set,
    # this needs go first as we need unmodified model. Also set to
    # preequilibration fixpars first as this is where initial states would be
    # computed otherwise
    model.setFixedParameters(edata.fixedParametersPreequilibration)
    edata.x0 = model.getInitialStates()
    edata.sx0 = model.getInitialStateSensitivities()

    # perturb model initial states
    model.setInitialStates(rdata['x_ss'] * 4)
    model.setInitialStateSensitivities(rdata['sx_ss'].flatten() / 2)

    # set ExpData plist
    edata.plist = model.getParameterList()
    # perturb model parameter list
    model.setParameterList([
        i for i in reversed(model.getParameterList())
    ])

    # set ExpData parameters
    edata.parameters = model.getParameters()
    # perturb model parameters
    model.setParameters(tuple(
        p * 2 for p in model.getParameters()
    ))

    # set ExpData pscale
    edata.pscale = model.getParameterScale()
    # perturb model pscale, needs to be done after getting parameters,
    # otherwise we will mess up parameter value
    model.setParameterScale(amici.parameterScalingFromIntVector([
        amici.ParameterScaling_log10
        if scaling == amici.ParameterScaling_none
        else amici.ParameterScaling_none
        for scaling in model.getParameterScale()
    ]))

    rdata_edata = amici.runAmiciSimulation(
        model, solver, edata
    )
    for variable in ['x', 'sx']:
        assert np.isclose(
            rdata[variable][0, :],
            rdata_edata[variable][0, :],
            1e-6, 1e-6
        ).all(), variable
