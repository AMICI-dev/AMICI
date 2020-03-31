import amici
import pytest
import os


@pytest.fixture(scope="session")
def sbml_model_presimulation_module():
    """Fixture for model_presimulation import"""

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
