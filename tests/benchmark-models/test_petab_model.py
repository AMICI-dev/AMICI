#!/usr/bin/env python3

"""
Simulate a PEtab problem and compare results to reference values
"""

import argparse
import copy
import importlib
import logging
import os
import sys
from typing import List, Union

import numpy as np
import pandas as pd
import petab
import yaml
import amici
from amici.logging import get_logger
from amici.petab_import import import_petab_problem
from amici.petab_objective import (simulate_petab, rdatas_to_measurement_df,
                                   LLH, SLLH, RDATAS)
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
    parser.add_argument('-s', '--sllh', dest='sllh',
                        action='store_true',
                        help=(
                            'Check sensitivities from '
                            '`amici.petab_objective.simulate_petab` with '
                            'finite differences.'
                        ))
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


def check_sllh(
    petab_problem: petab.Problem,
    amici_model: amici.Model,
    epsilon: Union[float, List[float]],
    abs_error_tolerance: float = 1e-1,
    rel_error_tolerance: float = 1e-1,
) -> bool:
    """Check sensitivities from `amici.petab_objective.simulate_petab`.

    Args:
        petab_problem:
            The PEtab problem.
        amici_model:
            The AMICI model for the PEtab problem.
        epsilon:
            The finite difference.

    Returns:
        Whether the finite difference check passes.
    """
    epsilons = epsilon
    if not isinstance(epsilon, list):
        epsilons = [epsilon]

    problem_parameters = dict(zip(
        petab_problem.x_ids,
        petab_problem.x_nominal,
    ))

    amici_solver = amici_model.getSolver()
    amici_solver.setSensitivityOrder(amici.SensitivityOrder_first)
    #amici_solver.setAbsoluteTolerance(1e-8)
    #amici_solver.setRelativeTolerance(1e-8)

    def simulate(
        vector,
        petab_problem=petab_problem,
        amici_model=amici_model,
        amici_solver=amici_solver,
    ):
        return amici.petab_objective.simulate_petab(
            petab_problem=petab_problem,
            amici_model=amici_model,
            solver=amici_solver,
            problem_parameters=vector,
        )

    result = simulate(problem_parameters)
    llh = result[LLH]
    x_ids = list(result[SLLH])

    central_differences = {x_id: [] for x_id in x_ids}
    for x_id in x_ids:
        for epsilon in epsilons:
            backward_parameters = copy.deepcopy(problem_parameters)
            forward_parameters = copy.deepcopy(problem_parameters)

            backward_parameters[x_id] -= epsilon
            forward_parameters[x_id] += epsilon

            # FIXME remove? reduces errors but maybe some parameter can be
            # negative in some model, e.g. on linear scale
            if backward_parameters[x_id] < 0:
                continue

            backward_llh = simulate(backward_parameters)[LLH]
            forward_llh = simulate(forward_parameters)[LLH]

            fd = (forward_llh - backward_llh) / (2 * epsilon)
            central_differences[x_id].append(fd)

    min_fd = {
        x_id: np.nanmin(differences)
        for x_id, differences in central_differences.items()
    }

    fd = pd.Series(min_fd)
    sllh = pd.Series(result[SLLH])

    abs_error = (sllh - fd).abs()
    rel_error = (abs_error / (fd + min(epsilons)*1e-3)).abs()

    print(rel_error)
    breakpoint()
    if (
        not (abs_error < abs_error_tolerance).all()
        or
        not (rel_error < rel_error_tolerance).all()
    ):
        raise ValueError('Sensitivities appear incorrect.')


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
    #if args.model_directory:
    #    sys.path.insert(0, args.model_directory)
    #model_module = importlib.import_module(args.model_name)
    #amici_model = model_module.getModel()
    amici_model = import_petab_problem(
        petab_problem=problem,
        model_output_dir=args.model_directory or None,
        model_name=args.model_name,
    )

    res = simulate_petab(
        petab_problem=problem, amici_model=amici_model,
        log_level=logging.DEBUG)
    rdatas = res[RDATAS]
    llh = res[LLH]

    # create simulation PEtab table
    sim_df = rdatas_to_measurement_df(rdatas=rdatas, model=amici_model,
                                      measurement_df=problem.measurement_df)
    sim_df.rename(columns={petab.MEASUREMENT: petab.SIMULATION}, inplace=True)

    if args.sllh:
        check_sllh(
            petab_problem=problem,
            amici_model=amici_model,
            epsilon=[1e-1, 1e-3, 1e-5, 1e-7, 1e-9],
        )

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
