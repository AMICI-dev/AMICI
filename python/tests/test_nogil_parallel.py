import os
import functools
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import amici
from amici.sbml_import import SbmlImporter
from amici.testing import TemporaryDirectoryWinSafe as TemporaryDirectory
import pytest

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


def _simulate(_: int, model_dir: str) -> int:
    model_module = amici.import_model_module("nogil_model", model_dir)
    model = model_module.getModel()
    model.setTimepoints(np.linspace(0, 2e5, 200001))
    solver = model.getSolver()
    rdata = amici.runAmiciSimulation(model, solver)
    return int(rdata.status)


def test_parallel_simulation_threading():
    assert os.cpu_count() >= 2, "requires at least two CPU cores"

    with TemporaryDirectory() as outdir:
        sbml_importer = SbmlImporter(SBML_EXAMPLE)
        sbml_importer.sbml2amici(model_name="nogil_model", output_dir=outdir)

        with ThreadPoolExecutor(max_workers=2) as pool:
            simulate = functools.partial(_simulate, model_dir=outdir)
            statuses = list(pool.map(simulate, range(100)))

    assert statuses == [amici.AMICI_SUCCESS] * 100
