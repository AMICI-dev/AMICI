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
        - {   eta      ,    t >= delta
    """
    # Model components
    model_name = 'piecewise'
    species = ['x_1', 'x_2']
    initial_assignments = {
        'x_1': 'zeta',
        'x_2': 'epsilon',
    }
    rate_rules = {
        'x_1': (
            'piecewise( alpha * x_1, time <    x_2, 0) + '
            'piecewise(- beta * x_1, time >= delta, 0)'
        ),
        'x_2': (
            'piecewise( gamma * x_2, time <  delta, 0) + '
            'piecewise(   eta      , time >= delta, 0)'
        ),
    }
    parameters = {
        'alpha': float(np.log(2)),
        'beta': float(np.log(4)),
        'gamma': float(np.log(3)),
        'delta': 1,
        'eta': 0.5,
        'zeta': 0.25,
        'epsilon': 100,
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
    def x_1(t, alpha, beta, gamma, delta, eta, zeta, epsilon):
        event_time = (
            (np.exp(gamma * delta) + delta * eta) /
                          (1 - eta)                      # noqa
        )
        if t < event_time:
            return zeta * np.exp(alpha * t)
        else:
            return zeta * np.exp(alpha * event_time - beta*(t - event_time))

    def x_2(t, alpha, beta, gamma, delta, eta, zeta, epsilon):
        event_time = delta
        if t < event_time:
            # TODO confirm usage of epsilon here
            return epsilon * np.exp(gamma*t)
        else:
            # TODO confirm usage of epsilon here
            return epsilon * np.exp(gamma*delta) + eta*(t-delta)
    result_expected = np.array([
        [x_1(t, **parameters) for t in timepoints],
        [x_2(t, **parameters) for t in timepoints],
    ])

    model.setTimepoints(timepoints)
    solver = model.getSolver()
    rdata = runAmiciSimulation(model, solver=solver)
    result_test = np.array([
        rdata['x'][:, 0],
        rdata['x'][:, 1],
    ])
    # The AMICI simulation matches the analytical solution.
    np.testing.assert_almost_equal(result_test, result_expected)
    # TODO test sensitivities directly
