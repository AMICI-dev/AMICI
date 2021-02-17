"""Pytest configuration for SBML test suite"""

import re
import sys
from typing import List

import pytest


# stores passed SBML semantic test suite IDs
passed_ids = []


def parse_selection(selection_str: str) -> List[int]:
    """
    Parse comma-separated list of integer ranges, return selected indices as
    integer list

    Valid input e.g.: "1", "1,3", "-3,4,6-7"
    """
    indices = []
    for group in selection_str.split(','):
        if not re.match(r'^(?:-?\d+)|(?:\d+(?:-\d+))$', group):
            print("Invalid selection", group)
            sys.exit()
        spl = group.split('-')
        if len(spl) == 1:
            indices.append(int(spl[0]))
        elif len(spl) == 2:
            begin = int(spl[0]) if spl[0] else 0
            end = int(spl[1])
            indices.extend(range(begin, end + 1))
    return indices


def pytest_addoption(parser):
    """Add pytest CLI options"""
    parser.addoption("--cases", help="Test cases to run")


def pytest_generate_tests(metafunc):
    """Parameterize tests"""

    # Run for all SBML semantic test suite cases
    if "test_number" in metafunc.fixturenames:
        # Get CLI option
        cases = metafunc.config.getoption("cases")
        if cases:
            # Run selected tests
            test_numbers = set(parse_selection(cases))
        else:
            # Run all tests
            test_numbers = set(range(1, 1781))

        # We skip this test due to NaNs in the Jacobian
        test_numbers -= {1395}

        metafunc.parametrize("test_number", test_numbers)


def pytest_sessionfinish(session, exitstatus):
    """Process test results"""
    global passed_ids
    terminalreporter = session.config.pluginmanager.get_plugin(
        'terminalreporter')
    terminalreporter.ensure_newline()
    # parse test names to get passed case IDs (don't know any better way to
    # access fixture values)
    from testSBMLSuite import format_test_id
    passed_ids = [format_test_id(_) for _ in passed_ids]
    if passed_ids:
        write_passed_tags(passed_ids, terminalreporter)
    terminalreporter.ensure_newline()


def write_passed_tags(passed_ids, out=sys.stdout):
    """Write tags of passed SBML semantic test cases"""
    passed_tags = set()
    from testSBMLSuite import get_tags_for_test
    for test_id in passed_ids:
        passed_tags |= get_tags_for_test(test_id)

    out.write("At least one test with the following tags has passed:\n")
    out.write('  ' + '\n  '.join(passed_tags))


def pytest_runtest_logreport(report: "TestReport") -> None:
    """Collect test case IDs of passed SBML semantic test suite cases"""
    if report.when == 'call'\
            and report.outcome == 'passed'\
            and '::test_sbml_testsuite_case[' in report.nodeid:
        test_case_id = re.sub(r'^.*::test_sbml_testsuite_case\[(\d+)].*$',
                              r'\1', report.nodeid)
        passed_ids.append(test_case_id)
