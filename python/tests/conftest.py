import copy
import importlib
import os
import sys

import amici
import pytest
from amici.pysb_import import pysb2amici

pysb = pytest.importorskip("pysb")


@pytest.fixture(scope="session")
def sbml_example_presimulation_module():
    """SBML example_presimulation model module fixture"""

    sbml_file = os.path.join(os.path.dirname(__file__), '..',
                             'examples', 'example_presimulation',
                             'model_presimulation.xml')

    sbml_importer = amici.SbmlImporter(sbml_file)

    constant_parameters = ['DRUG_0', 'KIN_0']

    observables = amici.assignmentRules2observables(
        sbml_importer.sbml,  # the libsbml model object
        filter_function=lambda variable: variable.getName() == 'pPROT_obs'
    )
    outdir = 'test_model_presimulation'
    module_name = 'test_model_presimulation'
    sbml_importer.sbml2amici(
        model_name=module_name,
        output_dir=outdir,
        verbose=False,
        observables=observables,
        constant_parameters=constant_parameters)

    model_module = amici.import_model_module(module_name=module_name,
                                             module_path=outdir)
    return model_module


@pytest.fixture(scope="session")
def pysb_example_presimulation_module():
    """PySB example_presimulation model module fixture"""

    constant_parameters = ['DRUG_0', 'KIN_0']

    pysb.SelfExporter.cleanup()  # reset pysb
    pysb.SelfExporter.do_export = True

    model_path = os.path.join(os.path.dirname(__file__), '..',
                              'examples', 'example_presimulation')

    with amici.add_path(model_path):
        if 'createModelPresimulation' in sys.modules:
            importlib.reload(sys.modules['createModelPresimulation'])
            model_module = sys.modules['createModelPresimulation']
        else:
            model_module = importlib.import_module('createModelPresimulation')
    model = copy.deepcopy(model_module.model)

    model.name = 'test_model_presimulation_pysb'
    outdir_pysb = model.name
    pysb2amici(model, outdir_pysb, verbose=False,
               observables=['pPROT_obs'],
               constant_parameters=constant_parameters)

    with amici.add_path(outdir_pysb):
        model_module_pysb = importlib.import_module(outdir_pysb)

    return model_module_pysb
