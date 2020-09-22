"""Test functions that require special treatment in amici.sbml_import.

In particular, some of the functions tested here require an installation of
boost.
"""

import os
import shutil

import amici
import numpy as np
from scipy.special import loggamma
import pytest
from amici.gradient_check import check_derivatives


def assert_fun(x):
    assert x


@pytest.fixture
def model_special_likelihoods():
    """Test model for special likelihood functions."""
    # load sbml model
    sbml_file = os.path.join(os.path.dirname(__file__), '..',
                             'examples', 'example_steadystate',
                             'model_steadystate_scaled.xml')
    sbml_importer = amici.SbmlImporter(sbml_file)

    # define observables
    observables = {
        'o1': {'formula': '100*10^x1'},
        'o2': {'formula': '100*10^x1'},
    }

    # define different noise models
    noise_distributions = {
        'o1': 'binomial', 'o2': 'negative-binomial',
    }

    module_name = 'test_special_likelihoods'
    outdir = 'test_special_likelihoods'
    sbml_importer.sbml2amici(
        model_name=module_name,
        output_dir=outdir,
        observables=observables,
        constant_parameters=['k0'],
        noise_distributions=noise_distributions,
    )

    yield amici.import_model_module(module_name=module_name,
                                    module_path=outdir)

    shutil.rmtree(outdir, ignore_errors=True)


def test_special_likelihoods(model_special_likelihoods):
    """Test special likelihood functions."""
    model = model_special_likelihoods.getModel()
    model.setTimepoints(np.linspace(0, 60, 10))
    solver = model.getSolver()
    solver.setSensitivityOrder(amici.SensitivityOrder.first)

    # Test in region with positive density

    # run model once to create an edata
    rdata = amici.runAmiciSimulation(model, solver)
    edata = amici.ExpData(rdata, 0.001, 0)

    # make sure measurements are smaller for non-degenerate probability
    y = edata.getObservedData()
    y = tuple([val * np.random.uniform(0, 1) for val in y])
    edata.setObservedData(y)

    # set sigmas
    sigma = 0.2
    sigmas = sigma * np.ones(len(y))
    edata.setObservedDataStdDev(sigmas)

    # and now run for real and also compute likelihood values
    rdata = amici.runAmiciSimulations(model, solver, [edata])[0]

    # check if the values make overall sense
    assert np.isfinite(rdata['llh'])
    assert np.all(np.isfinite(rdata['sllh']))
    assert np.any(rdata['sllh'])

    rdata_df = amici.getSimulationObservablesAsDataFrame(
        model, edata, rdata, by_id=True)
    edata_df = amici.getDataObservablesAsDataFrame(
        model, edata, by_id=True)

    # check correct likelihood value
    llh_exp = - sum([
        binomial_nllh(edata_df['o1'], rdata_df['o1'], sigma),
        negative_binomial_nllh(edata_df['o2'], rdata_df['o2'], sigma),
    ])
    assert np.isclose(rdata['llh'], llh_exp)

    # check gradient
    for sensi_method in [amici.SensitivityMethod.forward,
                         amici.SensitivityMethod.adjoint]:
        solver = model.getSolver()
        solver.setSensitivityMethod(sensi_method)
        solver.setSensitivityOrder(amici.SensitivityOrder.first)
        check_derivatives(
            model, solver, edata, assert_fun, atol=1e-1, rtol=1e-1,
            check_least_squares=False)

    # Test for m > y, i.e. in region with 0 density

    rdata = amici.runAmiciSimulation(model, solver)
    edata = amici.ExpData(rdata, 0.001, 0)

    # make sure measurements are smaller for non-degenerate probability
    y = edata.getObservedData()
    y = tuple([val * np.random.uniform(0.5, 3) for val in y])
    edata.setObservedData(y)
    edata.setObservedDataStdDev(sigmas)

    # and now run for real and also compute likelihood values
    rdata = amici.runAmiciSimulations(model, solver, [edata])[0]

    # m > y -> outside binomial domain -> 0 density
    assert rdata['llh'] == -np.inf
    # check for non-informative gradient
    assert all(np.isnan(rdata['sllh']))


def binomial_nllh(m: np.ndarray, y: np.ndarray, p: float):
    if any(m > y):
        return np.inf
    return sum(- loggamma(y+1) + loggamma(m+1) + loggamma(y-m+1) \
               - m * np.log(p) - (y-m) * np.log(1-p))


def negative_binomial_nllh(m: np.ndarray, y: np.ndarray, p: float):
    r = y * (1-p) / p
    return sum(- loggamma(m+r) + loggamma(m+1) + loggamma(r)
               - r * np.log(1-p) - m * np.log(p))
