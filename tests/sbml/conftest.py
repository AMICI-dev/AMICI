"""Pytest configuration for SBML test suite"""

import json
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

# the independently-reported checks per SBML semantic test suite case
SIMULATION_CHECK = "test_sbml_testsuite_case"
SENSITIVITY_FORWARD_CHECK = "test_sbml_testsuite_case_sensitivity_forward"
SENSITIVITY_ADJOINT_CHECK = "test_sbml_testsuite_case_sensitivity_adjoint"
SENSITIVITY_CONSISTENCY_CHECK = (
    "test_sbml_testsuite_case_sensitivity_consistency"
)
CHECKS = (
    SIMULATION_CHECK,
    SENSITIVITY_FORWARD_CHECK,
    SENSITIVITY_ADJOINT_CHECK,
    SENSITIVITY_CONSISTENCY_CHECK,
)
# short suffix per check, used for the `results.json` field names
_CHECK_SUFFIXES = {
    SIMULATION_CHECK: "simulation",
    SENSITIVITY_FORWARD_CHECK: "sensitivity_forward",
    SENSITIVITY_ADJOINT_CHECK: "sensitivity_adjoint",
    SENSITIVITY_CONSISTENCY_CHECK: "sensitivity_consistency",
}

# stores passed SBML semantic test suite IDs, by check
passed_ids: dict[str, list[str]] = {check: [] for check in CHECKS}
# test tags we encountered (from the simulation check only -- that's what
# the SBML test suite's own tag-support semantics are about)
encountered_tags: set[str] = set()
# failed/skipped tests with error message, by check
failed_or_skipped_ids: dict[str, dict[str, str]] = {
    check: {} for check in CHECKS
}

SBML_SEMANTIC_CASES_DIR = (
    Path(__file__).parent / "sbml-test-suite" / "cases" / "semantic"
)

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
    global passed_ids
    terminalreporter = session.config.pluginmanager.get_plugin(
        "terminalreporter"
    )
    terminalreporter.ensure_newline()
    # parse test names to get passed case IDs (don't know any better way to
    # access fixture values)
    passed_ids = {
        check: [format_test_id(_) for _ in ids]
        for check, ids in passed_ids.items()
    }
    if any(passed_ids.values()) or any(failed_or_skipped_ids.values()):
        write_passed_tags(passed_ids[SIMULATION_CHECK], terminalreporter)
    terminalreporter.ensure_newline()


def write_passed_tags(passed_simulation_ids, out=sys.stdout):
    """Write tags of passed SBML semantic test cases

    Tag coverage (what the SBML test suite's own result database tracks)
    only concerns basic simulation support, not the separate sensitivity
    check -- so tags are derived from `passed_simulation_ids` alone.
    """
    passed_component_tags = set()
    passed_test_tags = set()

    for test_id in passed_simulation_ids:
        cur_component_tags, cur_test_tags = get_tags_for_test(test_id)
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

    result = {
        "supported_tags": sorted(passed_test_tags | passed_component_tags),
        "encountered_tags": sorted(encountered_tags),
    }
    for check in CHECKS:
        suffix = _CHECK_SUFFIXES[check]
        ids = (
            passed_simulation_ids
            if check == SIMULATION_CHECK
            else passed_ids[check]
        )
        result[f"passed_tests_{suffix}"] = sorted(ids)
        result[f"failed_or_skipped_{suffix}"] = {
            k: failed_or_skipped_ids[check][k]
            for k in sorted(failed_or_skipped_ids[check])
        }

    with open(RESULT_PATH / "results.json", "w") as f:
        json.dump(result, f, indent=2)


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
    check, test_case_id = match.group(1), match.group(2)

    if check == SIMULATION_CHECK:
        global encountered_tags

        component_tags, test_tags = get_tags_for_test(test_case_id)
        encountered_tags |= component_tags
        encountered_tags |= test_tags

    if report.outcome == "passed":
        passed_ids[check].append(test_case_id)
    else:
        failed_or_skipped_ids[check][test_case_id] = report.longreprtext


def get_tags_for_test(test_id: str) -> tuple[set[str], set[str]]:
    """Get sbml test suite tags for the given test ID

    Returns:
        Tuple of set of strings for componentTags and testTags
    """
    current_test_path = SBML_SEMANTIC_CASES_DIR / test_id
    info_file = current_test_path / f"{test_id}-model.m"
    with open(info_file) as f:
        component_tags = set()
        test_tags = set()
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


def format_test_id(test_id) -> str:
    """Format numeric to 0-padded string"""
    return f"{test_id:0>5}"
