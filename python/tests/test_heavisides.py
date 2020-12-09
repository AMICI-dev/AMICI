"""Tests for SBML events, including piecewise expressions."""
import importlib
import libsbml
import numpy as np
from pathlib import Path
import pytest
import sys

from amici import (
    AmiciModel,
    runAmiciSimulation,
    SbmlImporter,
    SensitivityMethod,
    SensitivityOrder
)

sbml_test_models = Path('sbml_test_models')
sbml_test_models_output_dir = sbml_test_models / 'amici_models'
sbml_test_models_output_dir.mkdir(parents=True, exist_ok=True)


@pytest.fixture(params=[
    'state_and_parameter_dependent_heavisides',
    'piecewise_with_boolean_operations',
    'piecewise_many_conditions',
])
def model(request):
    """Returns the requested AMICI model and analytical expressions."""
    (
        initial_assignments,
        parameters,
        rate_rules,
        species,
        timepoints,
        x_pected,
        sx_pected
    ) = get_model_definition(request.param)

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
    amici_model = create_amici_model(
        sbml_model=sbml_model,
        model_name=request.param,
    )
    amici_model.setTimepoints(timepoints)

    return amici_model, parameters, timepoints, x_pected, sx_pected


def test_models(model):
    amici_model, parameters, timepoints, x_pected, sx_pected = model

    result_expected_x = np.array([
        x_pected(t, **parameters)
        for t in timepoints
    ])
    result_expected_sx = np.array([
        sx_pected(t, **parameters)
        for t in timepoints
    ])

    # --- Test the state trajectories without sensitivities -------------------
    # Does the AMICI simulation match the analytical solution?
    solver = amici_model.getSolver()
    rdata = runAmiciSimulation(amici_model, solver=solver)
    solver.setAbsoluteTolerance(1e-15)
    np.testing.assert_almost_equal(rdata['x'], result_expected_x, decimal=5)

    # Show that we can do arbitrary precision here (test 8 digits)
    solver = amici_model.getSolver()
    solver.setAbsoluteTolerance(1e-15)
    solver.setRelativeTolerance(1e-12)
    rdata = runAmiciSimulation(amici_model, solver=solver)
    np.testing.assert_almost_equal(rdata['x'], result_expected_x, decimal=8)

    # --- Test the state trajectories with sensitivities ----------------------

    # Does the AMICI simulation match the analytical solution?
    solver = amici_model.getSolver()
    solver.setSensitivityOrder(SensitivityOrder.first)
    solver.setSensitivityMethod(SensitivityMethod.forward)
    solver.setAbsoluteTolerance(1e-15)
    rdata = runAmiciSimulation(amici_model, solver=solver)

    np.testing.assert_almost_equal(rdata['x'], result_expected_x, decimal=5)
    np.testing.assert_almost_equal(rdata['sx'], result_expected_sx, decimal=5)

    # Show that we can do arbitrary precision here (test 8 digits)
    solver = amici_model.getSolver()
    solver.setSensitivityOrder(SensitivityOrder.first)
    solver.setSensitivityMethod(SensitivityMethod.forward)
    solver.setAbsoluteTolerance(1e-15)
    solver.setRelativeTolerance(1e-13)
    solver.setAbsoluteToleranceFSA(1e-15)
    solver.setRelativeToleranceFSA(1e-13)
    rdata = runAmiciSimulation(amici_model, solver=solver)
    np.testing.assert_almost_equal(rdata['x'], result_expected_x, decimal=8)
    np.testing.assert_almost_equal(rdata['sx'], result_expected_sx, decimal=8)


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


def get_model_definition(model_name):
    if model_name == 'state_and_parameter_dependent_heavisides':
        return model_definition_state_and_parameter_dependent_heavisides()
    elif model_name == 'piecewise_with_boolean_operations':
        return model_definition_piecewise_with_boolean_operations()
    elif model_name == 'piecewise_many_conditions':
        return model_definition_piecewise_many_conditions()
    else:
        raise NotImplementedError(
            f'Model with name {model_name} is not implemented.'
        )


def model_definition_state_and_parameter_dependent_heavisides():
    """Test model for state- and parameter-dependent heavisides.

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

    # Analytical solution
    def x_pected(t, alpha, beta, gamma, delta, eta, zeta):
        # get x_1
        tau_1 = (np.exp(gamma * delta) - delta * eta) / (1 - eta)
        if t < tau_1:
            x_1 = zeta * np.exp(alpha * t)
        else:
            x_1 = zeta * np.exp(alpha * tau_1 - beta*(t - tau_1))

        # get x_2
        tau_2 = delta
        if t < tau_2:
            x_2 = np.exp(gamma*t)
        else:
            x_2 = np.exp(gamma*delta) + eta*(t-delta)

        return (x_1, x_2)

    def sx_pected(t, alpha, beta, gamma, delta, eta, zeta):
        # get sx_1, w.r.t. parameters
        tau_1 = (np.exp(gamma * delta) - delta * eta) / (1 - eta)
        if t < tau_1:
            sx_1_alpha = zeta * t * np.exp(alpha * t)
            sx_1_beta = 0
            sx_1_gamma = 0
            sx_1_delta = 0
            sx_1_eta = 0
            sx_1_zeta = np.exp(alpha * t)
        else:
            # Never trust Wolfram Alpha...
            sx_1_alpha = (
                zeta * tau_1 * np.exp(alpha * tau_1 - beta*(t - tau_1))
            )
            sx_1_beta = (
                zeta * (tau_1 - t)
                * np.exp(alpha * tau_1 - beta*(t - tau_1))
            )
            sx_1_gamma = (
                zeta * (alpha + beta) * delta * np.exp(gamma * delta)
                / (1 - eta)
                * np.exp(alpha * tau_1 - beta*(t - tau_1))
            )
            sx_1_delta = (
                zeta * (alpha + beta)
                * np.exp(alpha * tau_1 - beta*(t - tau_1))
                * (gamma * np.exp(gamma * delta) - eta)
                / (1 - eta)
            )
            sx_1_eta = (
                zeta * (alpha + beta)
                * (-delta * (1-eta) + np.exp(gamma * delta) - delta * eta)
                / (1 - eta)**2
                * np.exp(alpha * tau_1 - beta*(t - tau_1))
            )
            sx_1_zeta = np.exp(alpha * tau_1 - beta*(t - tau_1))

        # get sx_2, w.r.t. parameters
        tau_2 = delta
        if t < tau_2:
            sx_2_alpha = 0
            sx_2_beta = 0
            sx_2_gamma = t * np.exp(gamma*t)
            sx_2_delta = 0
            sx_2_eta = 0
            sx_2_zeta = 0
        else:
            sx_2_alpha = 0
            sx_2_beta = 0
            sx_2_gamma = delta * np.exp(gamma*delta)
            sx_2_delta = gamma*np.exp(gamma*delta) - eta
            sx_2_eta = t - delta
            sx_2_zeta = 0

        sx_1 = (sx_1_alpha, sx_1_beta, sx_1_gamma,
                sx_1_delta, sx_1_eta, sx_1_zeta)
        sx_2 = (sx_2_alpha, sx_2_beta, sx_2_gamma,
                sx_2_delta, sx_2_eta, sx_2_zeta)

        return np.array((sx_1, sx_2)).transpose()

    return (
        initial_assignments,
        parameters,
        rate_rules,
        species,
        timepoints,
        x_pected,
        sx_pected
    )


def model_definition_piecewise_with_boolean_operations():
    # Model components
    species = ['x_1']
    initial_assignments = {'x_1': 'x_1_0'}
    rate_rules = {
        'x_1': (
            'piecewise('
                '1, '                                       # noqa
                    '(alpha <= time && time < beta) || '    # noqa
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

    # Analytical solution
    def x_pected(t, x_1_0, alpha, beta, gamma, delta):
        """Test model for boolean operations in a piecewise condition.

        ODEs
        ----
        d/dt x_1:
            - { 1, (alpha <= t and t < beta) or (gamma <= t and t < delta)
            - { 0, otherwise
        """
        if t < alpha:
            return (x_1_0,)
        elif alpha <= t < beta:
            return (x_1_0 + (t - alpha),)
        elif beta <= t < gamma:
            return (x_1_0 + (beta - alpha),)
        elif gamma <= t < delta:
            return (x_1_0 + (beta - alpha) + (t - gamma),)
        else:
            return (x_1_0 + (beta - alpha) + (delta - gamma), )

    def sx_pected(t, x_1_0, alpha, beta, gamma, delta):
        # x0 is very simple...
        sx_x0 = 1
        sx_alpha = 0
        sx_beta = 0
        sx_gamma = 0
        sx_delta = 0

        if t >= alpha:
            sx_alpha = -1
        if t >= beta:
            sx_beta = 1
        if t >= gamma:
            sx_gamma = -1
        if t >= delta:
            sx_delta = 1

        sx = (sx_alpha, sx_beta, sx_gamma, sx_delta, sx_x0)

        return np.array((sx,)).transpose()

    return (
        initial_assignments,
        parameters,
        rate_rules,
        species,
        timepoints,
        x_pected,
        sx_pected
    )


def model_definition_piecewise_many_conditions():
    # Model components
    species = ['x_1']
    initial_assignments = {'x_1': 'x_1_0'}
    t_final = 5

    pieces = 'piecewise('
    for t in range(t_final):
        if t > 0:
            pieces += ', '
        if t % 2 == 1:
            pieces += f'1, time < {t + 1}'
        else:
            pieces += f'0, time < {t + 1}'
    pieces += ', 0)'
    rate_rules = {'x_1': pieces, }

    parameters = {
        'x_1_0': 1,
    }
    timepoints = np.linspace(0, t_final, 100)

    # Analytical solution
    def x_pected(t, x_1_0):
        """Test model for piecewise functions with many pieces.

        ODEs
        ----
        d/dt x_1:
            - { 1,    floor(t) is odd
            - { 0,    otherwise
        """
        if np.floor(t) % 2 == 1:
            return (x_1_0 + (np.floor(t)-1)/2 + (t-np.floor(t)), )
        else:
            return (x_1_0 + np.floor(t)/2, )

    def sx_pected(t, x_1_0):
        return np.array([[1, ], ])

    return (
        initial_assignments,
        parameters,
        rate_rules,
        species,
        timepoints,
        x_pected,
        sx_pected
    )
