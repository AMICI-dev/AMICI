"""PYSB model tests"""

import importlib
import logging
import os
import platform
import shutil
import pytest

pysb = pytest.importorskip("pysb")

import amici
import numpy as np
import pysb.examples
import pytest
from amici.pysb_import import pysb2amici
from amici import ParameterScaling, parameterScalingFromIntVector
from pysb.simulator import ScipyOdeSimulator

from amici.gradient_check import check_derivatives


def test_compare_to_sbml_import(pysb_example_presimulation_module,
                                sbml_example_presimulation_module):
    # -------------- PYSB -----------------

    model_pysb = pysb_example_presimulation_module.getModel()

    edata = get_data(model_pysb)

    rdata_pysb = get_results(model_pysb, edata)

    # -------------- SBML -----------------

    model_sbml = sbml_example_presimulation_module.getModel()

    rdata_sbml = get_results(model_sbml, edata)

    # check if preequilibration fixed parameters are correctly applied:
    for rdata, model, importer in zip([rdata_sbml, rdata_pysb],
                                      [model_sbml, model_pysb],
                                      ['sbml', 'pysb']):
        # check equilibrium fixed parameters
        assert np.isclose(
            [sum(rdata["x_ss"][[1, 3]]), sum(rdata["x_ss"][[2, 4]])],
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

    skip_attrs = ['ptr', 'preeq_t', 'numsteps', 'preeq_numsteps',
                  'numrhsevals', 'numerrtestfails', 'order', 'J', 'xdot',
                  'preeq_wrms', 'preeq_cpu_time', 'cpu_time',
                  'cpu_timeB', 'w']

    for field in rdata_pysb:
        if field in skip_attrs:
            continue

        if rdata_pysb[field] is None:
            assert rdata_sbml[field] is None, field
        elif rdata_sbml[field] is None:
            assert rdata_pysb[field] is None, field
        elif np.isnan(rdata_sbml[field]).all():
            assert np.isnan(rdata_pysb[field]).all(), field
        elif np.isnan(rdata_pysb[field]).all():
            assert np.isnan(rdata_sbml[field]).all(), field
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
    pysb = pytest.importorskip("pysb")

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
                compute_conservation_laws=compute_conservation_laws,
                observables=list(pysb_model.observables.keys())
            )

            amici_model_module = amici.import_model_module(pysb_model.name,
                                                           outdir)

            model_pysb = amici_model_module.getModel()

            model_pysb.setTimepoints(tspan)

            solver = model_pysb.getSolver()
            solver.setMaxSteps(int(1e6))
            solver.setAbsoluteTolerance(atol)
            solver.setRelativeTolerance(rtol)
            rdata = amici.runAmiciSimulation(model_pysb, solver)

            # check agreement of species simulation

            assert np.isclose(rdata['x'],
                              pysb_simres.species, 1e-4, 1e-4).all()

            if example not in ['fricker_2010_apoptosis', 'fixed_initial',
                               'bngwiki_egfr_simple_deletemolecules']:
                if example in ['tyson_oscillator', 'bax_pore_sequential',
                               'bax_pore', 'kinase_cascade',
                               'bngwiki_egfr_simple',
                               'bngwiki_enzymatic_cycle_mm',
                               'bngwiki_simple']:
                    solver.setAbsoluteTolerance(1e-14)
                    solver.setRelativeTolerance(1e-14)
                    epsilon = 1e-4
                else:
                    solver.setAbsoluteTolerance(1e-10)
                    solver.setRelativeTolerance(1e-10)
                    epsilon = 1e-3
                model_pysb.setParameterScale(parameterScalingFromIntVector([
                    ParameterScaling.log10 if p > 0 else ParameterScaling.none
                    for p in model_pysb.getParameters()
                ]))
                check_derivatives(model_pysb, solver,
                                  epsilon=epsilon,
                                  rtol=1e-2,
                                  atol=1e-2,
                                  skip_zero_pars=True)

            shutil.rmtree(outdir, ignore_errors=True)


def get_data(model):
    solver = model.getSolver()
    model.setTimepoints(np.linspace(0, 60, 61))
    model.setSteadyStateSensitivityMode(
        amici.SteadyStateSensitivityMode.simulationFSA
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
        amici.SteadyStateSensitivityMode.simulationFSA
    )
    return amici.runAmiciSimulation(model, solver, edata)


def test_names_and_ids(pysb_example_presimulation_module):
    model_pysb = pysb_example_presimulation_module.getModel()
    expected = {
        'ExpressionIds': (
            '__s2',
            '__s1',
            '__s5',
            'pPROT',
            'tPROT',
            'initProt',
            'initDrug',
            'initKin',
            'pPROT_obs'),
        'FixedParameterIds': ('DRUG_0', 'KIN_0'),
        'FixedParameterNames': ('DRUG_0', 'KIN_0'),
        'ObservableIds': ('pPROT_obs',),
        'ObservableNames': ('pPROT_obs',),
        'ParameterIds': (
            'PROT_0',
            'kon_prot_drug',
            'koff_prot_drug',
            'kon_prot_kin',
            'kphospho_prot_kin',
            'kdephospho_prot'
        ),
        'StateIds': ('__s0', '__s1', '__s2', '__s3', '__s4', '__s5'),
        'StateNames': (
            "PROT(kin=None, drug=None, phospho='u')",
            'DRUG(bound=None)',
            'KIN(bound=None)',
            "DRUG(bound=1) % PROT(kin=None, drug=1, phospho='u')",
            "KIN(bound=1) % PROT(kin=1, drug=None, phospho='u')",
            "PROT(kin=None, drug=None, phospho='p')"
        ),
    }
    # Names and IDs are the same here
    expected['ExpressionNames'] = expected['ExpressionIds']
    expected['ParameterNames'] = expected['ParameterIds']

    for field_name, cur_expected in expected.items():
        actual = getattr(model_pysb, f'get{field_name}')()
        assert actual == cur_expected
