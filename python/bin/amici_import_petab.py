#!/usr/bin/env python3

"""
Import a model in the PEtab (https://github.com/ICB-DCM/PEtab/) format into
AMICI
"""

import amici
import petab
import os
import time
import argparse
from colorama import init as init_colorama
from colorama import Fore


def parse_cli_args():
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(
        description='Import PEtab-format model into AMICI.')

    # General options:
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        help='More verbose output')
    parser.add_argument('-o', '--output-dir', dest='model_output_dir',
                        help='Name of the model directory to create')

    # Call with set of files
    parser.add_argument('-s', '--sbml', dest='sbml_file_name',
                        help='SBML model filename')
    parser.add_argument('-m', '--measurements', dest='measurement_file_name',
                        help='Measurement table')
    parser.add_argument('-c', '--conditions', dest='condition_file_name',
                        help='Conditions table')
    parser.add_argument('-p', '--parameters', dest='parameter_file_name',
                        help='Parameter table')

    # or with model name, following default naming
    parser.add_argument('-n', '--model-name', dest='model_name',
                        help='Model name where all files are in the working '
                        'directory and follow PEtab naming convention. '
                        'Specifying -[smcp] will override defaults')
    args = parser.parse_args()

    if args.model_name:
        if not args.sbml_file_name:
            args.sbml_file_name = petab.get_default_sbml_file_name(
                args.model_name)
        if not args.measurement_file_name:
            args.measurement_file_name = \
                petab.get_default_measurement_file_name(args.model_name)
        if not args.condition_file_name:
            args.condition_file_name = petab.get_default_condition_file_name(
                args.model_name)
        if not args.parameter_file_name:
            args.parameter_file_name = petab.get_default_parameter_file_name(
                args.model_name)

    if not args.model_name and \
            (not args.sbml_file_name
             or not args.condition_file_name
             or not args.measurement_file_name):
        parser.error('When not specifying a model name, sbml, '
                     'condition and measurement file must be specified')

    return args


def show_model_info(sbml_model):
    """Print some model quantities"""

    print('Species: ', len(sbml_model.getListOfSpecies()))
    print('Global parameters: ', len(sbml_model.getListOfParameters()))
    print('Reactions: ', len(sbml_model.getListOfReactions()))


def get_fixed_parameters(condition_file_name, sbml_model,
                         constant_species_to_parameters=True):
    """Determine, set and return fixed model parameters

    Parameters specified in `condition_file_name` are turned into constants.
    Only global SBML parameters are considered. Local parameters are ignored.

    Species which are marked constant within the SBML model will be turned
    into constant parameters *within* the given `sbml_model`.
    """

    fixed_parameters_df = petab.get_condition_df(condition_file_name)
    print(f'Condition table: {fixed_parameters_df.shape}')

    # column names are model parameter names that should be made constant
    fixed_parameters = list(fixed_parameters_df.columns)
    try:
        fixed_parameters.remove('conditionName')
    except ValueError:
        pass
    # must be unique
    assert(len(fixed_parameters) == len(set(fixed_parameters)))

    if constant_species_to_parameters:
        # Turn species which are marked constant in the SBML model into
        # parameters
        constant_species = amici.sbml_import.constantSpeciesToParameters(
            sbml_model)
        print("Constant species converted to parameters",
              len(constant_species))
        print("Non-constant species", len(sbml_model.getListOfSpecies()))

    # ... and append them to the list of fixed_parameters
    for species in constant_species:
        if species not in fixed_parameters:
            fixed_parameters.append(species)

    # Ensure mentioned parameters exist in the model. remove additional ones
    for fixed_parameter in fixed_parameters[:]:
        # check global parameters
        if not sbml_model.getParameter(fixed_parameter):
            print(f"{Fore.YELLOW}Parameter '{fixed_parameter}' provided in "
                  "condition table but not present in model.")
            fixed_parameters.remove(fixed_parameter)

    return fixed_parameters


def ensure_amici_compatible():
    minimum_amici_version = (0, 8, 4)
    current_amici_version = \
        tuple([int(x) for x in amici.__version__.split('.')])
    if current_amici_version < minimum_amici_version:
        raise RuntimeError(f'Must use AMICI version {minimum_amici_version}'
                           ' or newer')


def import_model(sbml_file, condition_file, model_name=None,
                 model_output_dir=None, verbose=True):
    """Import AMICI model"""

    ensure_amici_compatible()

    if model_name is None:
        model_name = os.path.splitext(os.path.split(sbml_file)[-1])[0]

    if model_output_dir is None:
        model_output_dir = os.path.join(os.getcwd(), model_name)

    if verbose:
        print(f"{Fore.GREEN}Importing model '{sbml_file}' using fixed "
              f"parameters file '{condition_file}'")
        print(f"{Fore.GREEN}Model name is '{model_name}' Writing model code "
              f"to '{model_output_dir}'")

    sbml_importer = amici.SbmlImporter(sbml_file)
    sbml_model = sbml_importer.sbml

    show_model_info(sbml_model)

    observables = petab.get_observables(sbml_importer.sbml)

    sigmas = petab.get_sigmas(sbml_importer.sbml)

    if verbose:
        print('Observables', len(observables))
        print('Sigmas', len(sigmas))
    assert(len(sigmas) == len(observables))

    fixed_parameters = get_fixed_parameters(condition_file, sbml_model)

    if verbose:
        print("Overall fixed parameters", len(fixed_parameters))
        print("Non-constant global parameters",
              len(sbml_model.getListOfParameters()) - len(fixed_parameters))

    # Create Python module from SBML model
    start = time.time()
    sbml_importer.sbml2amici(model_name,
                             output_dir=model_output_dir,
                             observables=observables,
                             constantParameters=fixed_parameters,
                             sigmas=sigmas,
                             allow_reinit_fixpar_initcond=False,
                            )
    end = time. time()
    if verbose:
        print(f"{Fore.GREEN}Model imported successfully in "
              f"{round(end - start, 2)}s")


def main():
    args = parse_cli_args()

    """
    sbml_file = 'PEtab/CS_Signalling_ERBB_RAS_AKT_petab.xml'
    condition_file = 'PEtab/conditions_petab.tsv'
    measurement_file = 'PEtab/measurements_petab.tsv'
    parameter_file = 'PEtab/parameters_petab.tsv'
    """

    init_colorama(autoreset=True)

    # First check for valid PEtab
    pp = petab.Problem(args.sbml_file_name,
                       args.condition_file_name,
                       args.measurement_file_name,
                       args.parameter_file_name)
    petab.lint_problem(pp)

    import_model(args.sbml_file_name, args.condition_file_name,
                 model_output_dir=args.model_output_dir, verbose=True)


if __name__ == '__main__':
    main()
