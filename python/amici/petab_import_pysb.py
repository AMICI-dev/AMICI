"""
PySB-PEtab Import
-----------------
Import a model in the PySB-adapted :mod:`petab`
(https://github.com/PEtab-dev/PEtab) format into AMICI.
"""

import logging
from pathlib import Path
from typing import Optional, Union

import petab
import pysb
import sympy as sp
from petab.C import (CONDITION_NAME, NOISE_FORMULA, OBSERVABLE_FORMULA)
from petab.models.pysb_model import PySBModel

from . import petab_import
from .logging import get_logger, log_execution_time, set_log_level

logger = get_logger(__name__, logging.WARNING)


def _add_observation_model(
        pysb_model: pysb.Model,
        petab_problem: petab.Problem
):
    """Extend PySB model by observation model as defined in the PEtab
    observables table"""

    # add any required output parameters
    local_syms = {sp.Symbol.__str__(comp): comp for comp in
                  pysb_model.components if
                  isinstance(comp, sp.Symbol)}
    for formula in [*petab_problem.observable_df[OBSERVABLE_FORMULA],
                    *petab_problem.observable_df[NOISE_FORMULA]]:
        sym = sp.sympify(formula, locals=local_syms)
        for s in sym.free_symbols:
            if not isinstance(s, pysb.Component):
                p = pysb.Parameter(str(s), 1.0, _export=False)
                pysb_model.add_component(p)
                local_syms[sp.Symbol.__str__(p)] = p

    # add observables and sigmas to pysb model
    for (observable_id, observable_formula, noise_formula) \
            in zip(petab_problem.observable_df.index,
                   petab_problem.observable_df[OBSERVABLE_FORMULA],
                   petab_problem.observable_df[NOISE_FORMULA]):
        obs_symbol = sp.sympify(observable_formula, locals=local_syms)
        if observable_id in pysb_model.expressions.keys():
            obs_expr = pysb_model.expressions[observable_id]
        else:
            obs_expr = pysb.Expression(observable_id, obs_symbol,
                                       _export=False)
            pysb_model.add_component(obs_expr)
        local_syms[observable_id] = obs_expr

        sigma_id = f"{observable_id}_sigma"
        sigma_symbol = sp.sympify(
            noise_formula,
            locals=local_syms
        )
        sigma_expr = pysb.Expression(sigma_id, sigma_symbol, _export=False)
        pysb_model.add_component(sigma_expr)
        local_syms[sigma_id] = sigma_expr


@log_execution_time('Importing PEtab model', logger)
def import_model_pysb(
        petab_problem: petab.Problem,
        model_output_dir: Optional[Union[str, Path]] = None,
        verbose: Optional[Union[bool, int]] = True,
        model_name: Optional[str] = None,
        **kwargs
) -> None:
    """
    Create AMICI model from PySB-PEtab problem

    :param petab_problem:
        PySB PEtab problem

    :param model_output_dir:
        Directory to write the model code to. Will be created if doesn't
        exist. Defaults to current directory.

    :param verbose:
        Print/log extra information.

    :param model_name:
        Name of the generated model module

    :param kwargs:
        Additional keyword arguments to be passed to
        :meth:`amici.pysb_import.pysb2amici`.
    """
    set_log_level(logger, verbose)

    logger.info("Importing model ...")

    observable_table = petab_problem.observable_df
    if not isinstance(petab_problem.model, PySBModel):
        raise ValueError("Not a PySB model")
    pysb_model = petab_problem.model.model

    # TODO create copy?
    _add_observation_model(pysb_model, petab_problem)

    # For pysb, we only allow parameters in the condition table
    # those must be pysb model parameters (either natively, or output
    # parameters from measurement or condition table that have been added in
    # PysbPetabProblem)
    model_parameters = [p.name for p in pysb_model.parameters]
    for x in petab_problem.condition_df.columns:
        if x == CONDITION_NAME:
            continue

        if x not in model_parameters:
            # TODO
            raise NotImplementedError(
                "For PySB PEtab import, only model parameters, but no states "
                "or compartments are allowed in the condition table. "
                f"Offending column: {x}"
            )

    # TODO
    constant_parameters = petab_import.get_fixed_parameters(
        petab_problem)

    if observable_table is None:
        observables = None
        sigmas = None
        noise_distrs = None
    else:
        observables = [expr.name for expr in pysb_model.expressions
                       if expr.name in observable_table.index]

        sigmas = {obs_id: f"{obs_id}_sigma" for obs_id in observables}

        noise_distrs = petab_import.petab_noise_distributions_to_amici(
            observable_table)

    from amici.pysb_import import pysb2amici
    pysb2amici(model=pysb_model,
               output_dir=model_output_dir,
               model_name=model_name,
               verbose=True,
               observables=observables,
               sigmas=sigmas,
               constant_parameters=constant_parameters,
               noise_distributions=noise_distrs,
               **kwargs)
