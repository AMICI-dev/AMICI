import os
import functools
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from amici.sbml_import import SbmlImporter
from amici.testing import TemporaryDirectoryWinSafe as TemporaryDirectory
import pytest
import subprocess
import sys

pytestmark = pytest.mark.skipif(
    not os.environ.get("AMICI_NOGIL"),
    reason="only run in nogil workflow",
)


SBML_EXAMPLE = (
    Path(__file__).parents[2]
    / "doc"
    / "examples"
    / "getting_started"
    / "model_steadystate_scaled.xml"
)


def _simulate(_: int, model_dir: str):
    import amici

    model_module = amici.import_model_module("nogil_model", model_dir)
    model = model_module.getModel()
    model.setTimepoints(np.linspace(0, 2e5, 200001))
    solver = model.getSolver()
    return amici.runAmiciSimulation(model, solver)


def test_import_warns_without_python_gil(monkeypatch):
    monkeypatch.delenv("PYTHON_GIL", raising=False)

    proc = subprocess.run(
        [sys.executable, "-W", "default", "-c", "import amici"],
        capture_output=True,
        text=True,
        env=os.environ,
        check=True,
    )

    assert "RuntimeWarning: The global interpreter lock" in proc.stderr


def test_parallel_simulation_threading(monkeypatch):
    monkeypatch.setenv("PYTHON_GIL", "0")

    import amici

    assert os.cpu_count() >= 2, "requires at least two CPU cores"

    with TemporaryDirectory() as outdir:
        sbml_importer = SbmlImporter(SBML_EXAMPLE)
        sbml_importer.sbml2amici(model_name="nogil_model", output_dir=outdir)

        with ThreadPoolExecutor(max_workers=2) as pool:
            simulate = functools.partial(_simulate, model_dir=outdir)
            rdatas = list(pool.map(simulate, range(100)))

    statuses = [int(r.status) for r in rdatas]
    assert statuses == [amici.AMICI_SUCCESS] * 100

    first = rdatas[0]
    for rdata in rdatas[1:]:
        assert np.array_equal(rdata.x, first.x)
