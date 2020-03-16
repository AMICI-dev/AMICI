"""Pytest configuration for PEtab test suite"""

from typing import List
import re
import sys
import petabtests


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
    parser.addoption("--petab-cases", help="Test cases to run")


def pytest_generate_tests(metafunc):
    """Parameterize tests"""

    # Run for all PEtab test suite cases
    if "case" in metafunc.fixturenames:
        # Get CLI option
        cases = metafunc.config.getoption("--petab-cases")
        if cases:
            # Run selected tests
            test_numbers = parse_selection(cases)
        else:
            # Run all tests
            test_numbers = petabtests.CASES_LIST

        metafunc.parametrize("case", test_numbers)
