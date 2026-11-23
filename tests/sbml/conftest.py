"""Pytest configuration for SBML test suite"""

import os
import re
import sys
from pathlib import Path
from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
    from _pytest.reports import TestReport

# ensure that the script directory is in the Python path
script_dir = Path(__file__).parent.resolve()
if str(script_dir) not in sys.path:
    sys.path.insert(0, str(script_dir))

from suite_results import SuiteResults, read_tags  # noqa: E402

# test function -> check name in `results.json`, per test case
_CHECKS_BY_TEST = {
    "test_sbml_testsuite_case": "simulation",
    "test_sbml_testsuite_case_sensitivity_forward": "sensitivity_forward",
    "test_sbml_testsuite_case_sensitivity_adjoint": "sensitivity_adjoint",
    "test_sbml_testsuite_case_sensitivity_consistency": (
        "sensitivity_consistency"
    ),
}

# results of this run (`supported_tags` are set at the end)
results = SuiteResults()

# checkout of https://github.com/sbmlteam/sbml-test-suite/
SBML_TEST_SUITE_ROOT = Path(
    os.environ.get("SBML_TEST_SUITE_ROOT")
    or Path(__file__).parent / "sbml-test-suite"
)
SBML_SEMANTIC_CASES_DIR = SBML_TEST_SUITE_ROOT / "cases" / "semantic"

RESULT_PATH = Path(__file__).parent / "amici-semantic-results"


@pytest.fixture(scope="session")
def result_path() -> Path:
    return RESULT_PATH


@pytest.fixture(scope="session")
def sbml_semantic_cases_dir() -> Path:
    """directory with sbml semantic test cases"""
    return SBML_SEMANTIC_CASES_DIR


def parse_selection(selection_str: str, last: int) -> list[int]:
    """
    Parse comma-separated list of integer ranges, return selected indices as
    integer list

    Valid input e.g.: "1", "1,3", "-3,4,6-7"
    """
    indices = []
    for group in selection_str.split(","):
        if not re.match(r"^(?:-?\d+|\d+-\d*)$", group):
            print("Invalid selection", group)
            sys.exit()
        spl = group.split("-")
        if len(spl) == 1:
            indices.append(int(spl[0]))
        elif len(spl) == 2:
            begin = int(spl[0]) if spl[0] else 0
            end = int(spl[1]) if spl[1] else last
            indices.extend(range(begin, end + 1))
    return indices


def get_all_semantic_case_ids():
    """Get iterator over test sorted IDs of all cases in the SBML semantic
    suite"""
    pattern = re.compile(r"\d{5}")
    return sorted(
        str(x.name)
        for x in SBML_SEMANTIC_CASES_DIR.iterdir()
        if pattern.match(x.name)
    )


def pytest_addoption(parser):
    """Add pytest CLI options"""
    parser.addoption("--cases", help="Test cases to run")


def pytest_generate_tests(metafunc):
    """Parameterize tests"""

    # Run for all SBML semantic test suite cases
    if "test_id" in metafunc.fixturenames:
        # Get CLI option
        cases = metafunc.config.getoption("cases")
        if cases:
            # Run selected tests
            last_id = int(list(get_all_semantic_case_ids())[-1])
            test_numbers = sorted(set(parse_selection(cases, last_id)))
            test_ids = list(map(format_test_id, test_numbers))
        else:
            # Run all tests
            test_ids = get_all_semantic_case_ids()
        # `scope="session"` lets the session-scoped `compiled_case` fixture
        # (which imports/compiles the model -- expensive, must happen at
        # most once per test_id) depend on `test_id` at all; without it,
        # pytest raises `ScopeMismatch`. `xdist_group` keeps both this
        # test_id's simulation and sensitivity test nodes on the same
        # xdist worker -- required (not just an optimization) whenever
        # running with `-n`: `compiled_case`'s cache is per-worker-process,
        # so if the two nodes for one test_id land on different workers,
        # each independently recompiles the same model into the same
        # on-disk directory, which is both wasteful and racy. This only
        # takes effect when running with `--dist=loadgroup`.
        metafunc.parametrize(
            "test_id",
            [
                pytest.param(t, marks=pytest.mark.xdist_group(name=t))
                for t in test_ids
            ],
            scope="session",
        )


def pytest_sessionfinish(session, exitstatus):
    """Process test results"""
    terminalreporter = session.config.pluginmanager.get_plugin(
        "terminalreporter"
    )
    terminalreporter.ensure_newline()
    if any(results.passed.values()) or any(results.failed_or_skipped.values()):
        write_passed_tags(terminalreporter)
    terminalreporter.ensure_newline()


def write_passed_tags(out=sys.stdout):
    """Write tags of passed SBML semantic test cases

    Tag coverage (what the SBML test suite's own result database tracks)
    only concerns basic simulation support, not the separate sensitivity
    check -- so tags are derived from the simulation check alone.
    """
    passed_component_tags = set()
    passed_test_tags = set()

    for test_id in results.passed["simulation"]:
        cur_component_tags, cur_test_tags = read_tags(
            SBML_SEMANTIC_CASES_DIR / test_id, test_id
        )
        passed_component_tags |= cur_component_tags
        passed_test_tags |= cur_test_tags

    if passed_component_tags:
        out.write(
            "\nAt least one test with the following component tags has passed:\n"
        )
        out.write("  " + "\n  ".join(sorted(passed_component_tags)))

    if passed_test_tags:
        out.write(
            "\n\nAt least one test with the following test tags has passed:\n"
        )
        out.write("  " + "\n  ".join(sorted(passed_test_tags)))

    results.supported_tags = passed_test_tags | passed_component_tags
    results.to_json(RESULT_PATH / "results.json")


def pytest_runtest_logreport(report: "TestReport") -> None:
    """Collect test case IDs of passed SBML semantic test suite cases"""
    if report.when != "call":
        return
    match = re.search(
        r"::(test_sbml_testsuite_case"
        r"(?:_sensitivity_(?:forward|adjoint|consistency))?)\[(\d+)\]",
        report.nodeid,
    )
    if not match:
        return
    check = _CHECKS_BY_TEST[match.group(1)]
    test_case_id = format_test_id(match.group(2))

    if check == "simulation":
        component_tags, test_tags = read_tags(
            SBML_SEMANTIC_CASES_DIR / test_case_id, test_case_id
        )
        results.encountered_tags |= component_tags | test_tags

    if report.outcome == "passed":
        results.passed[check].add(test_case_id)
    else:
        results.failed_or_skipped[check][test_case_id] = report.longreprtext


def format_test_id(test_id) -> str:
    """Format numeric to 0-padded string"""
    return f"{test_id:0>5}"
