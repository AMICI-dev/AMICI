"""Tests related to amici.petab_import"""

import libsbml
import pandas as pd
import pytest
from amici.testing import skip_on_valgrind

petab = pytest.importorskip("petab.v1", reason="Missing petab")


@pytest.fixture
def simple_sbml_model():
    """Create a basic SBML model for test_get_fixed_parameters"""

    document = libsbml.SBMLDocument(3, 1)
    model = document.createModel()
    model.setId("simple_sbml_model")
    model.setTimeUnits("second")
    model.setExtentUnits("mole")
    model.setSubstanceUnits("mole")

    for par_idx in range(1, 6):
        p = model.createParameter()
        p.setId(f"p{par_idx}")
        p.setValue(par_idx)

    c = model.createCompartment()
    c.setId("c1")

    s = model.createSpecies()
    s.setId("x1")
    s.setConstant(True)
    s.setInitialConcentration(1.0)
    s.setCompartment(c.getId())

    return document, model


@pytest.fixture()
def get_fixed_parameters_model():
    """Create test SBML model for test_get_fixed_parameters"""
    ant_model = """
    p1 = 1
    p2 = 2
    p3 = 3
    p4 = 4
    p5 = 5
    p6 = 3^2
    p7 = p6
    p8 = 8
    p8' = 1
    p9 := p8
    """
    from amici.antimony_import import antimony2sbml

    sbml_str = antimony2sbml(ant_model)
    sbml_doc = libsbml.SBMLReader().readSBMLFromString(sbml_str)
    return sbml_doc, sbml_doc.getModel()


@skip_on_valgrind
def test_get_fixed_parameters(get_fixed_parameters_model):
    """Check for correct identification of fixed parameters:

    p1: fixed (via condition table)
    p2: (so far) not fixed (parametric override in condition table)
    p3: fixed (via parameter table `estimate=0`)
    p4: not fixed (via parameter table `estimate=1`)
    p5: fixed (implicitly, because not listed as estimated)
    p6: fixed (implicitly, because not listed as estimated
        initial assignment is a number)
    p7: not fixed (initial assignment is not a number)
    p8: not fixed (rate rule target)
    p9: not fixed (assignment rule target)
    """
    from amici.petab.sbml_import import (
        _get_fixed_parameters_sbml as get_fixed_parameters,
    )
    from petab.v1.models.sbml_model import SbmlModel

    sbml_doc, sbml_model = get_fixed_parameters_model
    condition_df = petab.get_condition_df(
        pd.DataFrame(
            {
                petab.CONDITION_ID: ["condition0"],
                "p1": [1.0],
                "p2": ["p1"],
            }
        )
    )
    parameter_df = petab.get_parameter_df(
        pd.DataFrame(
            {petab.PARAMETER_ID: ["p3", "p4"], petab.ESTIMATE: [0, 1]}
        )
    )
    print(condition_df)
    print(parameter_df)
    petab_problem = petab.Problem(
        model=SbmlModel(sbml_model),
        parameter_df=parameter_df,
        condition_df=condition_df,
    )
    assert set(get_fixed_parameters(petab_problem)) == {
        "p1",
        "p3",
        "p5",
        "p6",
    }

    assert set(
        get_fixed_parameters(
            petab_problem, non_estimated_parameters_as_constants=False
        )
    ) == {"p1", "p5", "p6"}


@skip_on_valgrind
def test_default_output_parameters(simple_sbml_model, tempdir):
    from amici.petab.petab_import import import_model
    from petab.v1.models.sbml_model import SbmlModel

    sbml_doc, sbml_model = simple_sbml_model
    condition_df = petab.get_condition_df(
        pd.DataFrame(
            {
                petab.CONDITION_ID: ["condition0"],
            }
        )
    )
    parameter_df = petab.get_parameter_df(
        pd.DataFrame({petab.PARAMETER_ID: [], petab.ESTIMATE: []})
    )
    observable_df = petab.get_observable_df(
        pd.DataFrame(
            {
                petab.OBSERVABLE_ID: ["obs1"],
                petab.OBSERVABLE_FORMULA: ["observableParameter1_obs1"],
                petab.NOISE_FORMULA: [1],
            }
        )
    )
    petab_problem = petab.Problem(
        model=SbmlModel(sbml_model),
        parameter_df=parameter_df,
        condition_df=condition_df,
        observable_df=observable_df,
    )

    sbml_importer = import_model(
        petab_problem=petab_problem,
        output_parameter_defaults={"observableParameter1_obs1": 1.0},
        compile=False,
        model_output_dir=tempdir,
    )
    assert (
        1.0
        == sbml_importer.sbml.getParameter(
            "observableParameter1_obs1"
        ).getValue()
    )

    with pytest.raises(ValueError):
        import_model(
            petab_problem=petab_problem,
            output_parameter_defaults={"nonExistentParameter": 1.0},
            compile=False,
            model_output_dir=tempdir,
        )
