import argparse

import petab.v1 as petab

from ..petab_import import import_model_sbml
from petab.v1.models.sbml_model import SbmlModel


def _parse_cli_args():
    """
    Parse command line arguments

    :return:
        Parsed CLI arguments from :mod:`argparse`.
    """
    parser = argparse.ArgumentParser(
        description="Import PEtab-format model into AMICI."
    )

    # General options:
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="More verbose output",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        dest="model_output_dir",
        help="Name of the model directory to create",
    )
    parser.add_argument(
        "--no-compile",
        action="store_false",
        dest="compile",
        help="Only generate model code, do not compile",
    )
    parser.add_argument(
        "--no-validate",
        action="store_false",
        dest="validate",
        help="Skip validation of PEtab files",
    )
    parser.add_argument(
        "--flatten",
        dest="flatten",
        default=False,
        action="store_true",
        help="Flatten measurement specific overrides of "
        "observable and noise parameters",
    )
    parser.add_argument(
        "--no-sensitivities",
        dest="generate_sensitivity_code",
        default=True,
        action="store_false",
        help="Skip generation of sensitivity code",
    )

    # Call with set of files
    group = parser.add_argument_group(
        "Providing individual PEtab tables *DEPRECATED*. "
        "Pass a PEtab yaml file instead."
    )

    group.add_argument(
        "-s", "--sbml", dest="sbml_file_name", help="SBML model filename"
    )
    group.add_argument(
        "-m",
        "--measurements",
        dest="measurement_file_name",
        help="Measurement table",
    )
    group.add_argument(
        "-c",
        "--conditions",
        dest="condition_file_name",
        help="Conditions table",
    )
    group.add_argument(
        "-p",
        "--parameters",
        dest="parameter_file_name",
        help="Parameter table",
    )
    group.add_argument(
        "-b",
        "--observables",
        dest="observable_file_name",
        help="Observable table",
    )

    parser.add_argument(
        "-y",
        "--yaml",
        dest="yaml_file_name_deprecated",
        help="PEtab YAML problem filename. *DEPRECATED* Pass the YAML file "
        "as positional argument instead.",
    )

    parser.add_argument(
        dest="yaml_file_name",
        help="PEtab YAML problem filename.",
        nargs="?",
    )

    parser.add_argument(
        "-n",
        "--model-name",
        dest="model_name",
        help="Name of the python module generated for the " "model",
    )

    args = parser.parse_args()
    if any(
        [
            args.sbml_file_name,
            args.condition_file_name,
            args.observable_file_name,
            args.measurement_file_name,
            args.parameter_file_name,
        ]
    ):
        print(
            "WARNING: Passing individual tables to amico_import_petab is "
            "deprecated, please pass a PEtab YAML file instead."
        )
    if (
        not args.yaml_file_name and not args.yaml_file_name_deprecated
    ) and not all(
        (
            args.sbml_file_name,
            args.condition_file_name,
            args.observable_file_name,
            args.measurement_file_name,
            args.parameter_file_name,
        )
    ):
        parser.error(
            "When not specifying a model name or YAML file, then "
            "SBML, condition, observable, measurement and parameter file must "
            "be specified."
        )

    if args.yaml_file_name_deprecated:
        print(
            "WARNING: -y/--yaml is deprecated. Pass the YAML file as "
            "positional argument instead."
        )
        args.yaml_file_name = args.yaml_file_name_deprecated

    return args


def _main():
    """
    Command line interface to import a model in the PEtab
    (https://github.com/PEtab-dev/PEtab/) format into AMICI.
    """
    args = _parse_cli_args()

    if args.yaml_file_name:
        pp = petab.Problem.from_yaml(args.yaml_file_name)
    else:
        pp = petab.Problem(
            model=SbmlModel.from_file(args.sbml_file_name),
            condition_df=petab.get_condition_df(args.condition_file_name),
            measurement_df=petab.get_measurement_df(
                args.measurement_file_name
            ),
            parameter_df=petab.get_parameter_df(args.parameter_file_name),
            observable_df=petab.get_observable_df(args.observable_file_name),
        )

    # Check for valid PEtab before potentially modifying it
    if args.validate:
        petab.lint_problem(pp)

    if args.flatten:
        petab.flatten_timepoint_specific_output_overrides(pp)

    import_model_sbml(
        model_name=args.model_name,
        petab_problem=pp,
        model_output_dir=args.model_output_dir,
        compile=args.compile,
        generate_sensitivity_code=args.generate_sensitivity_code,
        verbose=args.verbose,
        validate=False,
    )


if __name__ == "__main__":
    _main()
