#!/usr/bin/env python3
"""Generate HDF5 file with expected results for the C++ tests."""

import subprocess
from pathlib import Path

import amici
import h5py
from amici import (
    readModelDataFromHDF5,
    readSimulationExpData,
    readSolverSettingsFromHDF5,
    runAmiciSimulation,
    writeReturnData,
)
from amici.testing.models import (
    import_model_calvetti,
    import_model_dirac,
    import_model_events,
    import_model_jakstat,
    import_model_nested_events,
    import_model_neuron,
    import_model_robertson,
    import_model_steadystate,
)

repo_root = Path(__file__).parents[2]
outfile = Path(__file__).parent / "expected_results_py.h5"


def handle_model(id_: str, model: amici.Model):
    # write model settings to h5
    script = repo_root / "tests" / "generateTestConfig" / f"example_{id_}.py"
    subprocess.run([script, str(outfile)])

    # read test case IDs
    with h5py.File(outfile, "r") as f:
        cases = list(f[f"model_{id_}"])
    print(id_, cases)

    # generate expected results for each case
    for case in cases:
        # create a new model instance for each case, to ensure no interference
        model = model.module.get_model()
        readModelDataFromHDF5(
            str(outfile), model.get(), f"/model_{id_}/{case}/options"
        )
        solver = model.getSolver()
        readSolverSettingsFromHDF5(
            str(outfile), solver, f"/model_{id_}/{case}/options"
        )
        # read ExpData if data/ exists
        try:
            edata = readSimulationExpData(
                str(outfile), f"/model_{id_}/{case}/data", model.get()
            )
        except RuntimeError:
            edata = None

        rdata = runAmiciSimulation(model, solver, edata)
        writeReturnData(
            rdata._swigptr.get(), str(outfile), f"/model_{id_}/{case}/results"
        )


def main():
    """Generate expected results for the C++ tests."""
    handle_model("dirac", import_model_dirac())
    handle_model("events", import_model_events())
    handle_model("neuron", import_model_neuron())
    handle_model("calvetti", import_model_calvetti())
    handle_model("robertson", import_model_robertson())
    handle_model("jakstat_adjoint", import_model_jakstat())
    handle_model("nested_events", import_model_nested_events())
    handle_model("steadystate", import_model_steadystate())


if __name__ == "__main__":
    main()
