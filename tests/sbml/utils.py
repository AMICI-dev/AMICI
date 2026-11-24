from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import amici
import amici.sim.jax
import amici.sim.sundials as asd
import diffrax
import jax.numpy as jnp
import libsbml as sbml
import numpy as np
import optimistix
import pandas as pd
from amici.constants import SymbolId
from amici.importers.sbml import SbmlImporter
from amici.sim.jax import JAXModel
from amici.sim.jax.petab import (
    DEFAULT_CONTROLLER_SETTINGS,
    DEFAULT_ROOT_FINDER_SETTINGS,
)
from numpy.testing import assert_allclose
from sbml_cases import CaseSettings, find_model_file


def verify_results(
    settings: CaseSettings,
    rdata: Mapping[str, Any],
    expected: pd.DataFrame,
    importer: SbmlImporter,
    parameters: dict[str, float],
    expression_ids: Sequence[str],
) -> pd.DataFrame:
    """Verify simulation results against the expected results.

    :param rdata: The simulation results (`ts`, `y`, `w`).
    :param importer: The importer the simulated model was created with.
    :param parameters: Parameter ID -> value of the simulated model.
    :param expression_ids: Expression IDs of the simulated model.
    :return: The simulated values of `settings.variables` and time.
    """
    # collect states
    simulated = pd.DataFrame(
        rdata["y"],
        columns=[
            obs["name"]
            for obs in importer._symbols[SymbolId.OBSERVABLE].values()
        ],
    )
    simulated["time"] = rdata["ts"]

    # collect parameters
    parameter_data = {
        par_id: rdata["ts"] * 0 + value for par_id, value in parameters.items()
    }

    expression_data = {}

    for expr_idx, expr_id in enumerate(expression_ids):
        if expr_id.startswith("flux_"):
            new_key = expr_id.removeprefix("flux_")
        else:
            new_key = expr_id
            if expr_id in simulated.columns:
                continue  # skip if already present
        expression_data[new_key] = rdata["w"][:, expr_idx]

    # consolidated concatenation instead of columnwise insert to avoid data fragmentation in test 01395
    simulated = pd.concat(
        [
            simulated,
            pd.DataFrame(expression_data),
            pd.DataFrame(parameter_data),
        ],
        axis=1,
    )

    # SBML test suite case 01308 defines species with initialAmount and
    # hasOnlySubstanceUnits="true", but then request results as concentrations.
    # We only need to convert species that have only substance units
    concentration_species = [
        str(state_id)
        for state_id, state in {
            **importer._symbols[SymbolId.SPECIES],
            **importer._symbols[SymbolId.ALGEBRAIC_STATE],
        }.items()
        if str(state_id) in settings.concentrations
        and state.get("amount", False)
    ]
    amounts_to_concentrations(
        concentration_species, importer, simulated, settings.concentrations
    )

    concentrations_to_amounts(
        settings.amounts, importer, simulated, settings.concentrations
    )

    # simulated may contain `object` dtype columns and `expected` may
    # contain `np.int64` columns, so we cast everything to `np.float64`.
    for variable in settings.variables:
        expectation = expected[variable].astype(np.float64).values
        try:
            actual = simulated[variable].astype(np.float64).values
        except KeyError as e:
            raise KeyError(f"Missing simulated value for `{variable}`") from e
        assert_allclose(
            actual,
            expectation,
            settings.atol,
            settings.rtol,
            equal_nan=True,
            err_msg=f"Mismatch for {variable}",
        )

    return simulated[settings.variables + ["time"]]


def amounts_to_concentrations(
    amount_species: list[str],
    wrapper: SbmlImporter,
    simulated: pd.DataFrame,
    requested_concentrations: list[str],
) -> None:
    """
    Convert AMICI simulated amounts to concentrations

    Convert from concentration to amount:
    C=n/V
    n=CV (multiply by V)
    Convert from amount to concentration:
    n=CV
    C=n/V (divide by V)
    Dividing by V is equivalent to multiplying the reciprocal by V, then taking
    the reciprocal.
    This allows for the reuse of the concentrations_to_amounts method...
    """
    for species in amount_species:
        simulated.loc[:, species] = 1 / simulated.loc[:, species]
        concentrations_to_amounts(
            [species], wrapper, simulated, requested_concentrations
        )
        simulated.loc[:, species] = 1 / simulated.loc[:, species]


def concentrations_to_amounts(
    amount_species: list[str],
    wrapper: SbmlImporter,
    simulated: pd.DataFrame,
    requested_concentrations: list[str],
) -> None:
    """Convert AMICI simulated concentrations to amounts"""
    for species in amount_species:
        s = wrapper.sbml_model.getElementBySId(species)
        # Skip species that are marked to only have substance units since
        # they are already simulated as amounts
        if not isinstance(s, sbml.Species):
            continue

        is_amt = s.getHasOnlySubstanceUnits()
        comp = s.getCompartment()
        # Compartments and parameters that are treated as species do not
        # exist within a compartment.
        # Species with OnlySubstanceUnits don't have to be converted as long
        # as we don't request concentrations for them. Only applies when
        # called from amounts_to_concentrations.
        if (
            is_amt and species not in requested_concentrations
        ) or comp is None:
            continue

        simulated.loc[:, species] *= simulated.loc[:, comp]


def apply_settings(
    settings: CaseSettings, solver: asd.Solver, model: asd.Model, test_id: str
) -> None:
    """Apply model and solver settings as specified in the test case"""
    model.set_timepoints(settings.timepoints)
    solver.set_max_steps(int(1e6))
    solver.set_relative_tolerance(settings.rtol / 1e4)

    if test_id == "01148":
        solver.set_absolute_tolerance(settings.atol / 1e6)
    else:
        solver.set_absolute_tolerance(settings.atol / 1e4)


def compile_model_jax(
    case_dir: Path, test_id: str, model_dir: Path
) -> tuple[JAXModel, SbmlImporter]:
    """Import the given test model as JAX model.

    :return: The JAX model and the importer.
    """
    model_dir.mkdir(parents=True, exist_ok=True)
    importer = SbmlImporter(find_model_file(case_dir, test_id), jax=True)
    importer.sbml2jax(f"SBMLTest{test_id}_jax", output_dir=model_dir)
    model_module = amici.import_model_module(model_dir.name, model_dir.parent)
    return model_module.Model(), importer


def simulate_jax(
    model: JAXModel,
    parameters: jnp.ndarray,
    ts: np.ndarray | jnp.ndarray,
    atol: float,
    rtol: float,
    tol_factor: float = 1e2,
    root_finder: optimistix.AbstractRootFinder | None = None,
) -> tuple[jnp.ndarray, dict]:
    """Simulate a JAX model without data.

    :param parameters: The parameter values (`model.parameters` order).
    :param atol: The absolute tolerance for comparing against the expected
        results; the integration tolerance is `atol / tol_factor`.
    :param rtol: The relative tolerance, analogous to `atol`.
    :param root_finder: The root finder for events; defaults to
        `optimistix.Newton` with `DEFAULT_ROOT_FINDER_SETTINGS`.
    :return: The state trajectory and the simulation statistics, see
        `amici.sim.jax.JAXModel.simulate_condition`.
    """
    ts = jnp.asarray(ts, dtype=float)
    controller = diffrax.PIDController(
        rtol=rtol / tol_factor,
        atol=atol / tol_factor,
        pcoeff=DEFAULT_CONTROLLER_SETTINGS["pcoeff"],
        icoeff=DEFAULT_CONTROLLER_SETTINGS["icoeff"],
        dcoeff=DEFAULT_CONTROLLER_SETTINGS["dcoeff"],
    )
    if root_finder is None:
        root_finder = optimistix.Newton(**DEFAULT_ROOT_FINDER_SETTINGS)
    return model.simulate_condition(
        parameters,
        ts,
        jnp.array([]),
        jnp.zeros_like(ts),
        jnp.zeros_like(ts, dtype=int),
        jnp.zeros_like(ts, dtype=int),
        jnp.zeros((ts.shape[0], 0)),
        jnp.zeros((ts.shape[0], 0)),
        diffrax.Kvaerno5(),
        controller,
        root_finder,
        diffrax.DirectAdjoint(),
        diffrax.SteadyStateEvent(),
        2**10,
        ret=amici.sim.jax.ReturnValue.x,
    )
