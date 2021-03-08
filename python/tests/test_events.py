"""Tests for SBML events, including piecewise expressions."""
import importlib
import libsbml
import numpy as np
from pathlib import Path
import pytest
import sys
from scipy.linalg import expm
from copy import deepcopy


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
    'events_plus_heavisides', # 'nested_events'
])
def model(request):
    """Returns the requested AMICI model and analytical expressions."""
    (
        initial_assignments,
        parameters,
        rate_rules,
        species,
        events,
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
        events=events,
        to_file=sbml_test_models / ('event_test' + '.sbml'),
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
        sx_pected(t, parameters)
        for t in timepoints
    ])

    # --- Test the state trajectories without sensitivities -------------------
    # Does the AMICI simulation match the analytical solution?
    solver = amici_model.getSolver()
    solver.setAbsoluteTolerance(1e-15)
    rdata = runAmiciSimulation(amici_model, solver=solver)
    np.testing.assert_almost_equal(rdata['x'], result_expected_x, decimal=5)

    # Show that we can do arbitrary precision here (test 8 digits)
    solver = amici_model.getSolver()
    solver.setAbsoluteTolerance(1e-15)
    solver.setRelativeTolerance(1e-12)
    rdata = runAmiciSimulation(amici_model, solver=solver)
    np.testing.assert_almost_equal(rdata['x'], result_expected_x, decimal=8)

    # --- Test the state trajectories with sensitivities ----------------------
    # Show that we can do arbitrary precision here (test 8 digits)
    solver = amici_model.getSolver()
    solver.setSensitivityOrder(SensitivityOrder.first)
    solver.setSensitivityMethod(SensitivityMethod.forward)
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
        events,
        to_file: str = None,
):
    """Create an SBML model from simple definitions.

    See the model definitions and usage in :py:func:`model` for example input.

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

    for event_id, event_def in events.items():
        event = model.createEvent()
        event.setId(event_id)
        event.setName(event_id)
        event.setUseValuesFromTriggerTime(True)
        trigger = event.createTrigger()
        trigger.setMath(libsbml.parseL3Formula(event_def['trigger']))
        trigger.setPersistent(True)
        trigger.setInitialValue(True)
        assignment = event.createEventAssignment()
        assignment.setVariable(event_def['target'])
        assignment.setMath(libsbml.parseL3Formula(event_def['assignment']))

    if to_file:
        libsbml.writeSBMLToFile(
            document,
            str(to_file),
        )

    # Need to return document, else AMICI throws an error.
    # (possibly due to garbage collection?)
    return document, model


def get_model_definition(model_name):
    if model_name == 'events_plus_heavisides':
        return model_definition_events_plus_heavisides()
#    elif model_name == 'simple_event':
#        return model_definition_simple_event()
    else:
        raise NotImplementedError(
            f'Model with name {model_name} is not implemented.'
        )


def model_definition_events_plus_heavisides():
    """Test model for state- and parameter-dependent heavisides.

    ODEs
    ----
    d/dt x_1:
        - {            0,    t <  delta
        - { -alpha * x_1,    t >= delta
    d/dt x_2:
        - beta * x_1 - gamma * x_2
    d/dt x_3:
        - {     -eta * x_3,    t <  zeta
        - { -eta * x_3 + 1,    t >= zeta

    Events:
    -------
    event_1:
        trigger: k1 - x_3
        bolus: [[ -x_3 / 2],
                [        0],
                [        0]]
    event_2:
        trigger: t - zeta
        bolus: [[        0],
                [        0],
                [ zeta / 3]]
    """
    # Model components
    species = ['x_1', 'x_2', 'x_3']
    initial_assignments = {
        'x_1': 'k1',
        'x_2': 'k2',
        'x_3': 'k3',
    }
    rate_rules = {
        'x_1': 'piecewise( -alpha * x_1, time >= delta, 0)',
        'x_2': 'beta * x_1 - gamma * x_2',
        'x_3': '-eta * x_3 + piecewise( 1, time >= zeta, 0)',
    }
    parameters = {
        'k1':  2,
        'k2': 0.01,
        'k3': 5,
        'alpha': 2,
        'beta': 3,
        'gamma': 2,
        'delta': 3,
        'eta': 1,
        'zeta': 5,
    }
    events = {
        'event_1': {
            'trigger': 'x_3 < k1',
            'target': 'x_1',
            'assignment': 'x_1 - x_3 / 2'
        },
        'event_2': {
            'trigger': 'time >= zeta',
            'target': 'x_3',
            'assignment': 'x_3 + zeta / 3'
        }
    }
    timepoints = np.linspace(0, 8, 400)

    # Analytical solution
    def x_pected(t, k1, k2, k3, alpha, beta, gamma, delta, eta, zeta):
        # The system reads dx/dt = Ax + b
        # x0 = (k1, k2, k3)
        x0 = np.array([[k1], [k2], [k3]])

        # gather event time points
        event_1_time = (np.log(k3) - np.log(k1)) / eta  # k1 > x3
        event_2_time = delta
        event_3_time = zeta

        def get_early_x(t):
            # compute dynamics
            if t < event_1_time:
                # Define A
                A = np.array([[0, 0, 0],
                              [beta, -gamma, 0],
                              [0, 0, -eta]])
                exp_At = expm(t * A)
                return np.matmul(exp_At, x0)

            elif t <= event_2_time:
                # "simulate" until first event
                A = np.array([[0, 0, 0],
                              [beta, -gamma, 0],
                              [0, 0, -eta]])
                exp_At = expm(event_1_time * A)
                x1 = np.matmul(exp_At, x0)
                # apply bolus
                delta_x = np.array([[float(-x1[2, :] / 2)], [0], [0]])
                x1 += delta_x
                # "simulate" on
                exp_At = expm((t - event_1_time) * A)
                return np.matmul(exp_At, x1)

        if t < event_2_time:
            x = get_early_x(t).flatten()
        elif t < event_3_time:
            x2 = get_early_x(event_2_time)

            A = np.array([[-alpha, 0, 0],
                          [beta, -gamma, 0],
                          [0, 0, -eta]])
            exp_At = expm((t - event_2_time) * A)
            x = np.matmul(exp_At, x2).flatten()
        else:
            x2 = get_early_x(event_2_time)

            A = np.array([[-alpha, 0, 0],
                          [beta, -gamma, 0],
                          [0, 0, -eta]])
            exp_At = expm((event_3_time - event_2_time) * A)
            x3 = np.matmul(exp_At, x2)
            # apply bolus
            x3 += np.array([[0], [0], [zeta / 3]])

            hom_x = np.matmul(expm((t - event_3_time) * A), x3)
            inhom_x = [[0], [0],
                       [-np.exp(-eta * (t - event_3_time)) / (eta)
                        + 1 / (eta)]]

            x = (hom_x + inhom_x).flatten()

        return np.array(x)

    def sx_pected(t, parameters):
        # get sx, w.r.t. parameters, via finite differences
        sx = []

        for ip in parameters:
            eps = 1e-6
            perturbed_params = deepcopy(parameters)
            perturbed_params[ip] += eps
            sx_p = x_pected(t, **perturbed_params)
            perturbed_params[ip] -= 2*eps
            sx_m = x_pected(t, **perturbed_params)
            sx.append((sx_p - sx_m) / (2 * eps))

        return np.array(sx)

    return (
        initial_assignments,
        parameters,
        rate_rules,
        species,
        events,
        timepoints,
        x_pected,
        sx_pected
    )

