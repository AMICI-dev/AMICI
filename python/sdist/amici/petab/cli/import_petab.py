import argparse

import petab.v1 as petab

from ..petab_import import import_model_sbml


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
    if not args.yaml_file_name:
        parser.error("Missing PEtab yaml file.")

    return args


def _main():
    """
    Command line interface to import a model in the PEtab
    (https://github.com/PEtab-dev/PEtab/) format into AMICI.
    """
    args = _parse_cli_args()

    pp = petab.Problem.from_yaml(args.yaml_file_name)

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
