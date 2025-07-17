import os

import amici
import pytest
from test_pregenerated_models import (
    expected_results_matlab,
    expected_results_file,
    options_file,
    verify_simulation_results,
)
from amici.testing.models import import_model_neuron, import_model_events


@pytest.mark.skipif(
    os.environ.get("AMICI_SKIP_CMAKE_TESTS", "") == "TRUE",
    reason="skipping cmake based test",
)
@pytest.mark.parametrize(
    "model", [import_model_neuron(), import_model_events()]
)
def test_cases(model):
    expected_results = expected_results_matlab
    solver = model.getSolver()

    model_name = model.getName()
    # we need a different name for the model module to avoid collisions
    #  with the matlab-pregenerated models, but we need the old name for
    #  the expected results
    model_name = model_name.removesuffix("_py")

    for case in list(expected_results[model_name].keys()):
        if case.startswith("sensi2"):
            continue

        amici.readModelDataFromHDF5(
            options_file, model.get(), f"/{model_name}/{case}/options"
        )
        amici.readSolverSettingsFromHDF5(
            options_file, solver.get(), f"/{model_name}/{case}/options"
        )

        edata = None
        if "data" in expected_results[model_name][case].keys():
            edata = amici.readSimulationExpData(
                str(expected_results_file),
                f"/{model_name}/{case}/data",
                model.get(),
            )
        rdata = amici.runAmiciSimulation(model, solver, edata)

        verify_simulation_opts = dict()

        if model_name.startswith("model_neuron"):
            verify_simulation_opts["atol"] = 1e-5
            verify_simulation_opts["rtol"] = 1e-2

        verify_simulation_results(
            rdata,
            expected_results[model_name][case]["results"],
            **verify_simulation_opts,
        )
