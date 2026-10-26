"""Smoke test for an installed AMICI wheel.

Imports an SBML model, compiles it against the installed ``amici`` package
and simulates it with forward sensitivities. This checks that model
compilation works against the installed (relocated) package, not just that
``amici`` can be imported.

Usage::

    python test_wheel.py SBML_FILE

Run this from outside the AMICI source tree, so that the installed package
is imported.
"""

import sys
import tempfile
from pathlib import Path

import amici
import numpy as np
from amici.importers.sbml import SbmlImporter
from amici.sim.sundials import (
    SensitivityMethod,
    SensitivityOrder,
    run_simulation,
)


def main():
    sbml_file = Path(sys.argv[1]).resolve()

    amici_dir = Path(amici.__file__).parent
    print(f"Testing amici {amici.__version__} from {amici_dir}")
    if "site-packages" not in amici_dir.parts:
        raise AssertionError(
            f"amici was imported from {amici_dir}, not from an installed "
            "package."
        )

    model_name = "model_steadystate_scaled_wheel_test"
    # On Windows, the loaded model extension can't be deleted
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmpdir:
        model_dir = Path(tmpdir, model_name)
        SbmlImporter(sbml_file).sbml2amici(model_name, model_dir, verbose=True)
        model_module = amici.import_model_module(model_name, model_dir)

        model = model_module.get_model()
        model.set_timepoints(np.linspace(0, 60, 61))
        solver = model.create_solver()
        solver.set_sensitivity_order(SensitivityOrder.first)
        solver.set_sensitivity_method(SensitivityMethod.forward)
        rdata = run_simulation(model, solver)

        assert rdata.status == 0, rdata.status
        assert np.all(np.isfinite(rdata.x)), rdata.x
        assert np.all(np.isfinite(rdata.sx)), rdata.sx
        assert np.any(rdata.sx != 0.0)
        print(f"Simulation succeeded: x(t_end) = {rdata.x[-1]}")


if __name__ == "__main__":
    main()
