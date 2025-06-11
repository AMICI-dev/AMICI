"""AMICI-specific pytest fixtures.

This module is not intended to be imported directly, but to be used
in `pytest_plugins` in `conftest.py`.
"""

import pytest
from .testing import TemporaryDirectoryWinSafe


@pytest.fixture
def tempdir(*args, **kwargs):
    """Fixture that provides a temporary directory.

    The temporary directory will be removed after the test function
    completes. The fixture is scoped to the function, so each test function
    will get a new temporary directory. Therefore, this is only compatible with
    test cases and other function-scoped fixtures.

    Different from pytest's ``tmpdir`` fixture, this one will ignore certain
    cleanup errors on Windows, see :class:`TemporaryDirectoryWinSafe`.
    """
    with TemporaryDirectoryWinSafe(*args, **kwargs) as tmpdir:
        yield tmpdir
