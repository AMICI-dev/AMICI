import os

import pytest

pytestmark = pytest.mark.skipif(
    not os.environ.get("AMICI_NOGIL"),
    reason="only run in nogil workflow",
)


def test_import_warns_without_python_gil(monkeypatch):
    with pytest.raises(
        RuntimeWarning,
        match=r"The global interpreter lock \(GIL\) has been enabled to load module",
    ):
        import amici

        print(amici.__version__)
