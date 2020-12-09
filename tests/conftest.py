"""Pytest configuration for SBML test suite"""

from typing import List
import re
import sys

import _pytest.skipping as skipping


def show_skipped(tr, lines):
    for rep in tr.stats.get("skipped", []):
        pos = tr.config.cwd_relative_nodeid(rep.nodeid)
        reason = rep.longrepr[-1]
        if reason.startswith("Skipped: "):
            reason = reason[9:]
        verbose_word = skipping._get_report_str(tr.config, report=rep)
        lines.append(f"{verbose_word}\t{pos}: {reason}")


skipping.show_skipped = show_skipped


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
