"""Tests related to amici.petab_import"""

import libsbml
import pytest
import pandas as pd

petab = pytest.importorskip("petab")
amici_petab_import = pytest.importorskip("amici.petab_import")


@pytest.fixture
def simple_sbml_model():
    """Create a basic SBML model for test_get_fixed_parameters"""

    document = libsbml.SBMLDocument(3, 1)
    model = document.createModel()
    model.setTimeUnits("second")
    model.setExtentUnits("mole")
    model.setSubstanceUnits('mole')

    for par_idx in range(1, 6):
        p = model.createParameter()
        p.setId(f"p{par_idx}")
        p.setValue(par_idx)

    c = model.createCompartment()
    c.setId("c1")

    s = model.createSpecies()
    s.setId('x1')
    s.setConstant(True)
    s.setInitialConcentration(1.0)

    return document, model


def test_get_fixed_parameters(simple_sbml_model):
    """Check for correct identification of fixed parameters:

    p1: fixed (via condition table)
    p2: (so far) not fixed (parametric override in condition table)
    p3: fixed (via parameter table `estimate=0`)
    p4: not fixed (via parameter table `estimate=1`)
    p5: fixed (implicitly, because not listed as estimated)
    """
    sbml_doc, sbml_model = simple_sbml_model
    condition_df = petab.get_condition_df(
        pd.DataFrame({
            petab.CONDITION_ID: ["condition0"],
            "p1": [1.0],
            "p2": ["p1"],
        })
    )
    parameter_df = petab.get_parameter_df(
        pd.DataFrame({
            petab.PARAMETER_ID: ["p3", "p4"],
            petab.ESTIMATE: [0, 1]
        })
    )
    print(condition_df)
    print(parameter_df)
    petab_problem = petab.Problem(sbml_model=sbml_model,
                                  parameter_df=parameter_df,
                                  condition_df=condition_df)
    assert set(amici_petab_import.get_fixed_parameters(petab_problem)) \
        == {"p1", "p3", "p5"}

    # TODO: any non-estimated model parameter should be made constant
