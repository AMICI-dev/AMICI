#!/usr/bin/env python3

"""
Command line interface to import a model in the PEtab
(https://github.com/ICB-DCM/PEtab/) format into sAMICI
"""

import petab
import argparse

from amici.petab_import import import_model


def parse_cli_args():
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(
        description='Import PEtab-format model into AMICI.')

    # General options:
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        help='More verbose output')
    parser.add_argument('-o', '--output-dir', dest='model_output_dir',
                        help='Name of the model directory to create')
    parser.add_argument('--no-compile', action='store_false',
                        dest='compile',
                        help='Only generate model code, do not compile')

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


def main():
    args = parse_cli_args()

    # First check for valid PEtab
    pp = petab.Problem.from_files(
        sbml_file=args.sbml_file_name,
        condition_file=args.condition_file_name,
        measurement_file=args.measurement_file_name,
        parameter_file=args.parameter_file_name)
    petab.lint_problem(pp)

    import_model(model_name=args.model_name,
                 sbml_file=args.sbml_file_name,
                 condition_file=args.condition_file_name,
                 measurement_file=args.measurement_file_name,
                 model_output_dir=args.model_output_dir,
                 compile=args.compile,
                 verbose=True)


if __name__ == '__main__':
    main()
