"""Reading and writing SBML semantic test suite case files.

See https://github.com/sbmlteam/sbml-test-suite/blob/473e119dd57226c3a7a729d598f9007f06f781c3/cases/semantic/README.md.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


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


def read_tags(case_dir: Path, test_id: str) -> tuple[set[str], set[str]]:
    """Read the component tags and test tags of the given test case."""
    component_tags = set()
    test_tags = set()
    with open(case_dir / f"{test_id}-model.m") as f:
        for line in f:
            if line.startswith("testTags:"):
                test_tags = set(
                    re.split(r"[ ,:]", line[len("testTags:") :].strip())
                )
                test_tags.discard("")
            if line.startswith("componentTags:"):
                component_tags = set(
                    re.split(r"[ ,:]", line[len("componentTags:") :].strip())
                )
                component_tags.discard("")
            if test_tags and component_tags:
                return component_tags, test_tags
    print(f"No componentTags or testTags found for test case {test_id}.")
    return component_tags, test_tags


def write_result_file(
    simulated: pd.DataFrame, test_id: str, result_path: Path
) -> None:
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
