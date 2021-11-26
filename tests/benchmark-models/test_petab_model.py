#!/usr/bin/env python3

"""
Simulate a PEtab problem and compare results to reference values
"""

import argparse
import importlib
import logging
import os
import sys

import petab
import yaml
from amici.logging import get_logger
from amici.petab_objective import (simulate_petab, rdatas_to_measurement_df,
                                   LLH, RDATAS)
from petab.visualize import plot_petab_problem

logger = get_logger(f"amici.{__name__}", logging.WARNING)


def parse_cli_args():
    """Parse command line arguments

    Returns:
        Parsed CLI arguments from ``argparse``.
    """

    parser = argparse.ArgumentParser(
        description='Simulate PEtab-format model using AMICI.')

    # General options:
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        help='More verbose output')
    parser.add_argument('-c', '--check', dest='check', action='store_true',
                        help='Compare to reference value')
    parser.add_argument('-p', '--plot', dest='plot', action='store_true',
                        help='Plot measurement and simulation results')

    # PEtab problem
    parser.add_argument('-y', '--yaml', dest='yaml_file_name',
                        required=True,
                        help='PEtab YAML problem filename')

    # Corresponding AMICI model
    parser.add_argument('-m', '--model-name', dest='model_name',
                        help='Name of the AMICI module of the model to '
                             'simulate.', required=True)
    parser.add_argument('-d', '--model-dir', dest='model_directory',
                        help='Directory containing the AMICI module of the '
                             'model to simulate. Required if model is not '
                             'in python path.')

    parser.add_argument('-o', '--simulation-file', dest='simulation_file',
                        help='File to write simulation result to, in PEtab'
                        'measurement table format.')

    args = parser.parse_args()

    return args


def main():
    """Simulate the model specified on the command line"""

    args = parse_cli_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    logger.info(f"Simulating '{args.model_name}' "
                f"({args.model_directory}) using PEtab data from "
                f"{args.yaml_file_name}")

    # load PEtab files
    problem = petab.Problem.from_yaml(args.yaml_file_name)
    petab.flatten_timepoint_specific_output_overrides(problem)

    # load model
    if args.model_directory:
        sys.path.insert(0, args.model_directory)
    model_module = importlib.import_module(args.model_name)
    amici_model = model_module.getModel()
    amici_solver = amici_model.getSolver()

    if args.model_name == "Isensee_JCB2018":
        amici_solver.setAbsoluteTolerance(1e-12)
        amici_solver.setRelativeTolerance(1e-12)

    res = simulate_petab(
        petab_problem=problem, amici_model=amici_model,
        solver=amici_solver, log_level=logging.DEBUG)
    rdatas = res[RDATAS]
    llh = res[LLH]

    # create simulation PEtab table
    sim_df = rdatas_to_measurement_df(rdatas=rdatas, model=amici_model,
                                      measurement_df=problem.measurement_df)
    sim_df.rename(columns={petab.MEASUREMENT: petab.SIMULATION}, inplace=True)

    if args.simulation_file:
        sim_df.to_csv(index=False, sep="\t")

    if args.plot:
        try:
            # visualize fit
            axs = plot_petab_problem(petab_problem=problem, sim_data=sim_df)

            # save figure
            for plot_id, ax in axs.items():
                fig_path = os.path.join(args.model_directory,
                                        args.model_name + "_" + plot_id
                                        + "_vis.png")
                logger.info(f"Saving figure to {fig_path}")
                ax.get_figure().savefig(fig_path, dpi=150)

        except NotImplementedError:
            pass

    if args.check:
        references_yaml = os.path.join(os.path.dirname(__file__),
                                       "benchmark_models.yaml")
        with open(references_yaml) as f:
            refs = yaml.full_load(f)

        try:
            ref_llh = refs[args.model_name]["llh"]
            logger.info(f"Reference llh: {ref_llh}")

            if abs(ref_llh - llh) < 1e-3:
                logger.info(f"Computed llh {llh} matches reference "
                            f"{ref_llh}. Absolute difference is "
                            f"{ref_llh - llh}.")
            else:
                logger.error(f"Computed llh {llh} does not match reference "
                             f"{ref_llh}. Absolute difference is "
                             f"{ref_llh - llh}."
                             f" Relative difference is {llh / ref_llh}")
                sys.exit(1)
        except KeyError:
            logger.error("No reference likelihood found for "
                         f"{args.model_name} in {references_yaml}")


if __name__ == "__main__":
    main()
