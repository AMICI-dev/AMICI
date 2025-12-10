"""pytest configuration file"""

import copy
import importlib
import os
import sys
from pathlib import Path

import amici
import pytest
from amici import MeasurementChannel
from amici.testing import TemporaryDirectoryWinSafe as TemporaryDirectory

pytest_plugins = ["amici.testing.fixtures"]

REPO_ROOT = Path(__file__).parents[2]
EXAMPLES_DIR = REPO_ROOT / "doc" / "examples"
TEST_DIR = Path(__file__).parent
MODEL_STEADYSTATE_SCALED_XML = (
    EXAMPLES_DIR / "getting_started" / "model_steadystate_scaled.xml"
)
MODEL_PRESIMULATION_XML = (
    EXAMPLES_DIR / "example_presimulation" / "model_presimulation.xml"
)

os.environ.setdefault("AMICI_MODELS_ROOT", str(REPO_ROOT.absolute()))


@pytest.fixture(scope="session")
def sbml_example_presimulation_module():
    """SBML example_presimulation model module fixture"""

    sbml_importer = amici.SbmlImporter(MODEL_PRESIMULATION_XML)

    fixed_parameters = ["DRUG_0", "KIN_0"]

    observables = amici.assignment_rules_to_observables(
        sbml_importer.sbml_model,  # the libsbml model object
        filter_function=lambda variable: variable.getName() == "pPROT_obs",
    )
    module_name = "test_model_presimulation"

    with TemporaryDirectory(prefix=module_name) as outdir:
        sbml_importer.sbml2amici(
            model_name=module_name,
            output_dir=outdir,
            verbose=False,
            observation_model=observables,
            fixed_parameters=fixed_parameters,
        )

        yield amici.import_model_module(
            module_name=module_name, module_path=outdir
        )


@pytest.fixture(scope="session")
def pysb_example_presimulation_module():
    """PySB example_presimulation model module fixture"""
    pysb = pytest.importorskip("pysb")
    from amici.importers.pysb import pysb2amici

    fixed_parameters = ["DRUG_0", "KIN_0"]

    pysb.SelfExporter.cleanup()  # reset pysb
    pysb.SelfExporter.do_export = True

    model_path = MODEL_PRESIMULATION_XML.parent

    with amici.add_path(model_path):
        if "createModelPresimulation" in sys.modules:
            importlib.reload(sys.modules["createModelPresimulation"])
            model_module = sys.modules["createModelPresimulation"]
        else:
            model_module = importlib.import_module("createModelPresimulation")
    model = copy.deepcopy(model_module.model)

    model.name = "test_model_presimulation_pysb"

    with TemporaryDirectory(prefix=model.name) as outdir:
        pysb2amici(
            model,
            outdir,
            verbose=True,
            observation_model=[MeasurementChannel("pPROT_obs")],
            fixed_parameters=fixed_parameters,
        )

        yield amici.import_model_module(model.name, outdir)


@pytest.fixture(scope="session")
def model_units_module():
    sbml_file = TEST_DIR / "model_units.xml"
    module_name = "test_model_units"

    sbml_importer = amici.SbmlImporter(sbml_file)

    with TemporaryDirectory() as outdir:
        sbml_importer.sbml2amici(model_name=module_name, output_dir=outdir)

        yield amici.import_model_module(
            module_name=module_name, module_path=outdir
        )
