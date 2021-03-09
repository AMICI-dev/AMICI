"""Tests for SBML events, including piecewise expressions."""
import importlib
import libsbml
import numpy as np
from pathlib import Path
import sys

from amici import (
    AmiciModel,
    import_model_module,
    runAmiciSimulation,
    SbmlImporter,
    SensitivityMethod,
    SensitivityOrder
)


def create_amici_model(sbml_model, model_name) -> AmiciModel:
    """
    Import an sbml file and create an AMICI model from it
    """
    sbml_test_models = Path('sbml_test_models')
    sbml_test_models_output_dir = sbml_test_models / 'amici_models'
    sbml_test_models_output_dir.mkdir(parents=True, exist_ok=True)

    sbml_importer = SbmlImporter(sbml_model)
    output_dir = sbml_test_models_output_dir / model_name
    sbml_importer.sbml2amici(
        model_name=model_name,
        output_dir=str(output_dir)
    )

    model_module = import_model_module(model_name, str(output_dir.resolve()))
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


def check_trajectories_without_sensitivities(
        amici_model: AmiciModel,
        result_expected_x: np.ndarray,
):
    """
    Check whether the AMICI simulation matches a known solution
    (ideally an analytically calculated one).
    """

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


def check_trajectories_with_forward_sensitivities(
        amici_model: AmiciModel,
        result_expected_x: np.ndarray,
        result_expected_sx: np.ndarray,
):
    """
    Check whether the forward sensitivities of the AMICI simulation match
    a known solution (ideally an analytically calculated one).
    """

    # Show that we can do arbitrary precision here (test 8 digits)
    solver = amici_model.getSolver()
    solver.setAbsoluteTolerance(1e-15)
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
