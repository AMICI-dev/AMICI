from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import amici.sim.sundials as asd
import libsbml as sbml
import numpy as np
import pandas as pd
from amici.constants import SymbolId
from amici.importers.sbml import SbmlImporter
from numpy.testing import assert_allclose


def _parse_list(value: str) -> list[str]:
    """Parse a comma-separated settings value."""
    return [item for item in value.replace(" ", "").split(",") if item]


@dataclass
class CaseSettings:
    """Settings of an SBML semantic test suite case.

    See the `<case>-settings.txt` files.
    """

    start: float
    duration: float
    steps: int
    #: IDs of all variables for which results are expected
    variables: list[str]
    #: Absolute tolerance for comparing against the expected results
    atol: float
    #: Relative tolerance for comparing against the expected results
    rtol: float
    #: Species for which results are expected as amounts
    amounts: list[str]
    #: Species for which results are expected as concentrations
    concentrations: list[str]

    @classmethod
    def from_file(cls, case_dir: Path, test_id: str) -> CaseSettings:
        """Read the settings file of the given test case."""
        raw = {}
        with open(case_dir / f"{test_id}-settings.txt") as f:
            for line in f:
                if line.strip():
                    key, value = line.split(":")
                    raw[key] = value.strip()
        # start/duration/steps may be empty
        return cls(
            start=float(raw["start"] or 0),
            duration=float(raw["duration"] or 0),
            steps=int(raw["steps"] or 0),
            variables=_parse_list(raw["variables"]),
            atol=float(raw["absolute"]),
            rtol=float(raw["relative"]),
            amounts=_parse_list(raw["amount"]),
            concentrations=_parse_list(raw["concentration"]),
        )

    @property
    def timepoints(self) -> np.ndarray:
        """The output timepoints."""
        return np.linspace(
            self.start, self.start + self.duration, self.steps + 1
        )


def read_expected_results(case_dir: Path, test_id: str) -> pd.DataFrame:
    """Read the expected results of the given test case."""
    results = pd.read_csv(case_dir / f"{test_id}-results.csv", delimiter=",")
    return results.rename(
        columns={c: c.replace(" ", "") for c in results.columns},
    )


def verify_results(
    settings: CaseSettings,
    rdata,
    expected: pd.DataFrame,
    importer: SbmlImporter,
    parameters: dict[str, float],
    expression_ids: list[str],
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
    amount_species, wrapper, simulated, requested_concentrations
):
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
    amount_species, wrapper, simulated, requested_concentrations
):
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


def write_result_file(
    simulated: pd.DataFrame, test_id: str, result_path: Path
):
    """
    Create test result file for upload to
    http://raterule.caltech.edu/Facilities/Database

    Requires csv file with test ID in name and content of [time, Species, ...]
    """
    # TODO: only states are reported here, not compartments or parameters

    result_path.mkdir(parents=True, exist_ok=True)
    filename = result_path / f"{test_id}.csv"
    (
        simulated
        # see https://github.com/sbmlteam/sbml-test-suite/blob/473e119dd57226c3a7a729d598f9007f06f781c3/cases/semantic/README.md?plain=1#L211-L215
        .replace({np.nan: "NaN", np.inf: "INF", -np.inf: "-INF"}).to_csv(
            filename, index=False
        )
    )


def apply_settings(
    settings: CaseSettings, solver: asd.Solver, model: asd.Model, test_id: str
):
    """Apply model and solver settings as specified in the test case"""
    model.set_timepoints(settings.timepoints)
    solver.set_max_steps(int(1e6))
    solver.set_relative_tolerance(settings.rtol / 1e4)

    if test_id == "01148":
        solver.set_absolute_tolerance(settings.atol / 1e6)
    else:
        solver.set_absolute_tolerance(settings.atol / 1e4)


def find_model_file(current_test_path: Path, test_id: str) -> Path:
    """Find model file for the given test (guess filename extension)"""

    sbml_file = current_test_path / f"{test_id}-sbml-l3v2.xml"

    if not sbml_file.is_file():
        # fallback l3v1
        sbml_file = current_test_path / f"{test_id}-sbml-l3v1.xml"

    if not sbml_file.is_file():
        # fallback l2v5
        sbml_file = current_test_path / f"{test_id}-sbml-l2v5.xml"

    return sbml_file
