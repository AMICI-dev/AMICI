"""Results of SBML semantic test suite runs (`results.json`).

Only uses the standard library, so that results can be consolidated
without AMICI installed.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

#: The independently reported checks per test case
CHECKS = (
    "simulation",
    "sensitivity_forward",
    "sensitivity_adjoint",
    "sensitivity_consistency",
)


@dataclass
class SuiteResults:
    """Results of an SBML semantic test suite run."""

    #: Tags of all test cases run (simulation check only)
    encountered_tags: set[str] = field(default_factory=set)
    #: Tags of the test cases that passed the simulation check
    supported_tags: set[str] = field(default_factory=set)
    #: IDs of passed test cases, by check
    passed: dict[str, set[str]] = field(
        default_factory=lambda: {check: set() for check in CHECKS}
    )
    #: IDs of failed or skipped test cases -> message, by check
    failed_or_skipped: dict[str, dict[str, str]] = field(
        default_factory=lambda: {check: {} for check in CHECKS}
    )

    @property
    def unsupported_tags(self) -> set[str]:
        """Tags for which no test case passed the simulation check."""
        return self.encountered_tags - self.supported_tags

    def update(self, other: SuiteResults) -> None:
        """Add the results of another (partial) run."""
        self.encountered_tags |= other.encountered_tags
        self.supported_tags |= other.supported_tags
        for check in CHECKS:
            self.passed[check] |= other.passed[check]
            self.failed_or_skipped[check] |= other.failed_or_skipped[check]

    def to_json(self, path: Path) -> None:
        """Write the results to a JSON file."""
        result = {
            "supported_tags": sorted(self.supported_tags),
            "encountered_tags": sorted(self.encountered_tags),
        }
        for check in CHECKS:
            result[f"passed_tests_{check}"] = sorted(self.passed[check])
            result[f"failed_or_skipped_{check}"] = dict(
                sorted(self.failed_or_skipped[check].items())
            )
        with open(path, "w") as f:
            json.dump(result, f, indent=2)

    @classmethod
    def from_json(cls, path: Path) -> SuiteResults:
        """Read results from a JSON file written by `to_json`."""
        with open(path) as f:
            data = json.load(f)
        return cls(
            encountered_tags=set(data["encountered_tags"]),
            supported_tags=set(data["supported_tags"]),
            passed={
                check: set(data[f"passed_tests_{check}"]) for check in CHECKS
            },
            failed_or_skipped={
                check: data[f"failed_or_skipped_{check}"] for check in CHECKS
            },
        )
