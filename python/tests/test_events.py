"""Tests for SBML events, including piecewise expressions."""
import importlib
import libsbml
import numpy as np
from pathlib import Path
import sys

from amici import (
    AmiciModel,
    runAmiciSimulation,
    SbmlImporter,
)

sbml_test_models = Path('sbml_test_models')
sbml_test_models_output_dir = sbml_test_models / 'amici_models'
sbml_test_models_output_dir.mkdir(parents=True, exist_ok=True)


def create_amici_model(sbml_model, model_name) -> AmiciModel:
    sbml_importer = SbmlImporter(sbml_model)
    output_dir = sbml_test_models_output_dir / model_name
    sbml_importer.sbml2amici(
        model_name=model_name,
        output_dir=str(output_dir)
    )

    sys.path.insert(0, str(output_dir.resolve()))
    model_module = importlib.import_module(model_name)

    model = model_module.getModel()
    return model


def create_sbml_model(
        initial_assignments,
        parameters,
        rate_rules,
        species,
        to_file: str = None,
):
    """Create an SBML model from simple definitions.

    See usage in :py:func:`test_piecewise` for example input.

    The default initial concentration of species is `1.0`. This can currently
    be changed by specifying an initial assignment.
    """
    document = libsbml.SBMLDocument(3, 1)
    model = document.createModel()

    compartment = model.createCompartment()
    compartment.setId('compartment')
    compartment.setConstant(True)
    compartment.setSize(1)
    compartment.setSpatialDimensions(3)
    compartment.setUnits('dimensionless')

    for species_id in species:
        species = model.createSpecies()
        species.setId(species_id)
        species.setCompartment('compartment')
        species.setConstant(False)
        species.setSubstanceUnits('dimensionless')
        species.setBoundaryCondition(False)
        species.setHasOnlySubstanceUnits(False)
        species.setInitialConcentration(1.0)

    for target, formula in initial_assignments.items():
        initial_assignment = model.createInitialAssignment()
        initial_assignment.setSymbol(target)
        initial_assignment.setMath(libsbml.parseL3Formula(formula))

    for target, formula in rate_rules.items():
        rate_rule = model.createRateRule()
        rate_rule.setVariable(target)
        rate_rule.setMath(libsbml.parseL3Formula(formula))

    for parameter_id, parameter_value in parameters.items():
        parameter = model.createParameter()
        parameter.setId(parameter_id)
        parameter.setConstant(True)
        parameter.setValue(parameter_value)
        parameter.setUnits('dimensionless')

    if to_file:
        libsbml.writeSBMLToFile(
            document,
            str(to_file),
        )

    # Need to return document, else AMICI throws an error.
    # (possibly due to garbage collection?)
    return document, model


def test_piecewise():
    """Test model for piecewise functions in ODEs.

    ODEs
    ----
    d/dt x_1:
        - { alpha * x_1,    t <  x_2
        - { -beta * x_1,    t >= x_2
    d/dt x_2:
        - { gamma * x_2,    t <  delta
        - {         eta,    t >= delta
    """
    # Model components
    model_name = 'piecewise'
    species = ['x_1', 'x_2']
    initial_assignments = {
        'x_1': 'zeta',
    }
    rate_rules = {
        'x_1': 'piecewise( alpha * x_1, time < x_2,   -beta * x_1 )',
        'x_2': 'piecewise( gamma * x_2, time < delta,  eta        )',
    }
    parameters = {
        'alpha': float(np.log(2)),
        'beta': float(np.log(4)),
        'gamma': float(np.log(3)),
        'delta': 1,
        'eta': 0.5,
        'zeta': 0.25,
    }
    timepoints = np.linspace(0, 10, 100)

    # SBML model
    sbml_document, sbml_model = create_sbml_model(
        initial_assignments=initial_assignments,
        parameters=parameters,
        rate_rules=rate_rules,
        species=species,
        # uncomment `to_file` to save SBML model to file for inspection
        # to_file=sbml_test_models / (model_name + '.sbml'),
    )

    # AMICI model
    model = create_amici_model(
        sbml_model=sbml_model,
        model_name=model_name,
    )

    # Analytical solution
    def x_1(t, alpha, beta, gamma, delta, eta, zeta):
        event_time = (
            (np.exp(gamma * delta) - delta * eta) / (1 - eta)
        )
        if t < event_time:
            return zeta * np.exp(alpha * t)
        else:
            return zeta * np.exp(alpha * event_time - beta*(t - event_time))

    def x_2(t, alpha, beta, gamma, delta, eta, zeta):
        event_time = delta
        if t < event_time:
            return np.exp(gamma*t)
        else:
            return np.exp(gamma*delta) + eta*(t-delta)

    result_expected = np.array([
        [x_1(t, **parameters) for t in timepoints],
        [x_2(t, **parameters) for t in timepoints],
    ]).transpose()

    model.setTimepoints(timepoints)
    solver = model.getSolver()
    rdata = runAmiciSimulation(model, solver=solver)
    result_test = rdata['x']

    # The AMICI simulation matches the analytical solution.
    np.testing.assert_almost_equal(result_test, result_expected, decimal=5)
    # Show that we can do arbitrary precision here (test 8 digits)
    solver = model.getSolver()
    solver.setRelativeTolerance(1e-12)
    rdata = runAmiciSimulation(model, solver=solver)
    result_test = rdata['x']
    np.testing.assert_almost_equal(result_test, result_expected, decimal=8)

    # TODO test sensitivities directly


def test_piecewise_complex_condition():
    """Test model for piecewise functions in ODEs.

    ODEs
    ----
    d/dt x_1:
        - { 1,    (alpha <= t and t < beta) or (gamma <= t and t < delta)
        - { 0,    otherwise
    """
    # Model components
    model_name = 'piecewise_complex'
    species = ['x_1']
    initial_assignments = {'x_1': 'x_1_0'}
    rate_rules = {
        'x_1': (
            'piecewise('
                '1, '
                    '(alpha <= time && time < beta) || '
                    '(gamma <= time && time < delta), '
                '0'
            ')'
        ),
    }
    parameters = {
        'alpha': 1,
        'beta': 2,
        'gamma': 3,
        'delta': 4,
        'x_1_0': 1,
    }
    timepoints = np.linspace(0, 5, 100)

    # SBML model
    sbml_document, sbml_model = create_sbml_model(
        initial_assignments=initial_assignments,
        parameters=parameters,
        rate_rules=rate_rules,
        species=species,
        # uncomment `to_file` to save SBML model to file for inspection
        to_file=sbml_test_models / (model_name + '.sbml'),
    )

    # AMICI model
    model = create_amici_model(
        sbml_model=sbml_model,
        model_name=model_name,
    )

    # Analytical solution
    def x_1(t, x_1_0, alpha, beta, gamma, delta):
        if t < alpha:
            return x_1_0
        elif alpha <= t < beta:
            return x_1_0 + (t - alpha)
        elif beta <= t < gamma:
            return x_1_0 + (beta - alpha)
        elif gamma <= t < delta:
            return x_1_0 + (beta - alpha) + (t - gamma)
        else:
            return x_1_0 + (beta - alpha) + (delta - gamma)

    result_expected = np.array([
        [x_1(t, **parameters) for t in timepoints],
    ]).transpose()

    model.setTimepoints(timepoints)
    solver = model.getSolver()
    rdata = runAmiciSimulation(model, solver=solver)
    result_test = rdata['x']

    # The AMICI simulation matches the analytical solution.
    np.testing.assert_almost_equal(result_test, result_expected, decimal=5)
    # Show that we can do arbitrary precision here (test 8 digits)
    solver = model.getSolver()
    solver.setRelativeTolerance(1e-12)
    rdata = runAmiciSimulation(model, solver=solver)
    result_test = rdata['x']
    np.testing.assert_almost_equal(result_test, result_expected, decimal=8)

    # TODO test sensitivities directly
