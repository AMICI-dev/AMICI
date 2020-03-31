"PYSB model tests"

import copy
import importlib
import logging
import os
import platform
import sys

import amici
import numpy as np
import pysb.examples
import pytest
from amici.pysb_import import pysb2amici
from pysb.simulator import ScipyOdeSimulator


@pytest.fixture
def pysb_example_presimulation():
    """PySB example_presimulation model module fixture"""
    constant_parameters = ['DRUG_0', 'KIN_0']

    pysb.SelfExporter.cleanup()  # reset pysb
    pysb.SelfExporter.do_export = True

    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..',
                                    'examples',
                                    'example_presimulation'))
    if 'createModelPresimulation' in sys.modules:
        importlib.reload(sys.modules['createModelPresimulation'])
        model_module = sys.modules['createModelPresimulation']
    else:
        model_module = importlib.import_module('createModelPresimulation')
    model = copy.deepcopy(model_module.model)
    model.name = 'test_model_presimulation_pysb'
    outdir_pysb = model.name
    pysb2amici(model, outdir_pysb, verbose=False,
               observables=['pPROT_obs'],
               constant_parameters=constant_parameters)
    sys.path.insert(0, outdir_pysb)
    modelModulePYSB = importlib.import_module(outdir_pysb)

    return modelModulePYSB


def test_compare_to_sbml_import(pysb_example_presimulation,
                                sbml_model_presimulation_module):
    # -------------- PYSB -----------------

    model_pysb = pysb_example_presimulation.getModel()

    edata = get_data(model_pysb)

    rdata_pysb = get_results(model_pysb, edata)

    # -------------- SBML -----------------

    modelModuleSBML = sbml_model_presimulation_module
    model_sbml = modelModuleSBML.getModel()

    rdata_sbml = get_results(model_sbml, edata)

    # check if preequilibration fixed parameters are correctly applied:
    for rdata, model, importer in zip([rdata_sbml, rdata_pysb],
                                      [model_sbml, model_pysb],
                                      ['sbml', 'pysb']):
        # check equilibrium fixed parameters
        assert np.isclose(
            [
                sum(rdata["x_ss"][[1, 3]]),
                sum(rdata["x_ss"][[2, 4]])
            ],
            edata.fixedParametersPreequilibration,
            atol=1e-6, rtol=1e-6
        ).all(), f'{importer} preequilibration'

        # check equilibrium initial parameters
        assert np.isclose(
            sum(rdata["x_ss"][[0, 3, 4, 5]]),
            model.getParameterByName('PROT_0'),
            atol=1e-6, rtol=1e-6
        ), f'{importer} preequilibration'

        # check reinitialization with fixed parameter after
        # presimulation
        assert np.isclose(
            [rdata["x0"][1], rdata["x0"][2]],
            edata.fixedParameters,
            atol=1e-6, rtol=1e-6
        ).all(), f'{importer} presimulation'

    skip_attrs = ['ptr', 't_steadystate', 'numsteps', 'newton_numsteps',
                  'numrhsevals', 'numerrtestfails', 'order', 'J', 'xdot',
                  'wrms_steadystate', 'newton_cpu_time', 'cpu_time',
                  'cpu_timeB', 'w']

    for field in rdata_pysb:
        if field in skip_attrs:
            continue

        if rdata_pysb[field] is None:
            assert rdata_sbml[field] is None, field
        elif rdata_sbml[field] is None:
            assert rdata_pysb[field] is None, field
        else:
            assert np.isclose(
                rdata_sbml[field], rdata_pysb[field],
                atol=1e-6, rtol=1e-6
            ).all(), field


pysb_models = [
    'tyson_oscillator', 'robertson', 'expression_observables',
    'bax_pore_sequential', 'bax_pore', 'bngwiki_egfr_simple',
    'bngwiki_enzymatic_cycle_mm', 'bngwiki_simple', 'earm_1_0',
    'earm_1_3', 'move_connected', 'michment', 'kinase_cascade',
    'hello_pysb', 'fricker_2010_apoptosis', 'explicit',
    'fixed_initial',
]
custom_models = [
    'bngwiki_egfr_simple_deletemolecules',
]


@pytest.mark.parametrize('example', pysb_models + custom_models)
def test_compare_to_pysb_simulation(example):
    atol = 1e-8
    rtol = 1e-8

    with amici.add_path(os.path.dirname(pysb.examples.__file__)):
        with amici.add_path(os.path.join(os.path.dirname(__file__), '..',
                                         'tests', 'pysb_test_models')):

            if example == 'earm_1_3' \
                    and platform.sys.version_info[0] == 3 \
                    and platform.sys.version_info[1] < 7:
                return

            # load example

            pysb.SelfExporter.cleanup()  # reset pysb
            pysb.SelfExporter.do_export = True

            module = importlib.import_module(example)
            pysb_model = module.model
            pysb_model.name = pysb_model.name.replace('pysb.examples.', '')
            # avoid naming clash for custom pysb models
            pysb_model.name += '_amici'

            # pysb part

            tspan = np.linspace(0, 100, 101)
            sim = ScipyOdeSimulator(
                pysb_model,
                tspan=tspan,
                integrator_options={'rtol': rtol, 'atol': atol},
                compiler='python'
            )
            pysb_simres = sim.run()

            # amici part

            outdir = pysb_model.name

            if pysb_model.name in ['move_connected_amici']:
                with pytest.raises(Exception):
                    pysb2amici(pysb_model, outdir, verbose=logging.INFO,
                               compute_conservation_laws=True)
                compute_conservation_laws = False
            else:
                compute_conservation_laws = True

            pysb2amici(
                pysb_model,
                outdir,
                verbose=logging.INFO,
                compute_conservation_laws=compute_conservation_laws
            )
            sys.path.insert(0, outdir)

            amici_model_module = importlib.import_module(pysb_model.name)

            model_pysb = amici_model_module.getModel()

            model_pysb.setTimepoints(tspan)

            solver = model_pysb.getSolver()
            solver.setMaxSteps(int(1e5))
            solver.setAbsoluteTolerance(atol)
            solver.setRelativeTolerance(rtol)
            rdata = amici.runAmiciSimulation(model_pysb, solver)

            # check agreement of species simulation

            assert np.isclose(rdata['x'],
                              pysb_simres.species, 1e-4, 1e-4).all()


def get_data(model):
    solver = model.getSolver()
    model.setTimepoints(np.linspace(0, 60, 61))
    model.setSteadyStateSensitivityMode(
        amici.SteadyStateSensitivityMode_simulationFSA
    )

    rdata = amici.runAmiciSimulation(model, solver)
    edata = amici.ExpData(rdata, 0.1, 0.0)
    edata.t_presim = 2
    edata.fixedParameters = [10, 2]
    edata.fixedParametersPresimulation = [3, 2]
    edata.fixedParametersPreequilibration = [3, 0]
    edata.reinitializeFixedParameterInitialStates = True
    return edata


def get_results(model, edata):
    solver = model.getSolver()
    solver.setSensitivityOrder(1)
    edata.reinitializeFixedParameterInitialStates = True
    model.setTimepoints(np.linspace(0, 60, 61))
    model.setSteadyStateSensitivityMode(
        amici.SteadyStateSensitivityMode_simulationFSA
    )
    return amici.runAmiciSimulation(model, solver, edata)
