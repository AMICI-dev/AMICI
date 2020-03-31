"""Tests related to amici.petab_import"""

import libsbml
import pytest

amici_petab_import = pytest.importorskip("amici.petab_import")


@pytest.fixture
def simple_sbml_model():
    """Create a basic SBML model for test_constant_species_to_parameters"""

    document = libsbml.SBMLDocument(3, 1)
    model = document.createModel()
    model.setTimeUnits("second")
    model.setExtentUnits("mole")
    model.setSubstanceUnits('mole')

    s = model.createSpecies()
    s.setId('x1')
    s.setConstant(True)
    s.setInitialConcentration(1.0)

    # Add reaction to ensure species is replaced in reactions
    r = model.createReaction()
    r.setId('r1')
    # Add multiple instances to ensure all are removed
    for (coeff, name) in [(1, 'x1'), (1, 'x1')]:
        species_ref = r.createReactant()
        species_ref.setSpecies(name)
        species_ref.setStoichiometry(coeff)
    for (coeff, name) in [(1, 'x1'), (1, 'x1')]:
        species_ref = r.createProduct()
        species_ref.setSpecies(name)
        species_ref.setStoichiometry(coeff)
    for name in ['x1', 'x1']:
        species_ref = r.createModifier()
        species_ref.setSpecies(name)

    return document, model


def test_constant_species_to_parameters(simple_sbml_model):
    """Test conversion from species to constant parameters"""

    document, model = simple_sbml_model

    amici_petab_import.constant_species_to_parameters(model)

    assert len(list(model.getListOfParameters())) == 1
    assert len(list(model.getListOfSpecies())) == 0

    r = model.getReaction(0)
    assert len(list(r.getListOfReactants())) == 0
    assert len(list(r.getListOfProducts())) == 0
    assert len(list(r.getListOfModifiers())) == 0
