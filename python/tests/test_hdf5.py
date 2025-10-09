"""AMICI HDF5 I/O tests"""

import os
import random

import amici
import pytest


def _modify_solver_attrs(solver):
    # change to non-default values
    for attr in dir(solver):
        if not attr.startswith("set"):
            continue

        val = getattr(solver, attr.replace("set", "get"))()

        if isinstance(val, bool):
            cval = not val
        elif attr == "set_stability_limit_flag":
            cval = 0
        elif attr == "set_return_data_reporting_mode":
            cval = amici.RDataReporting.likelihood
        elif attr == "set_max_time":
            # default value is the maximum, must not add to that
            cval = random.random()
        elif attr == "set_constraints":
            cval = [1.0, 1.0]
        elif isinstance(val, int):
            cval = val + 1
        else:
            cval = val + random.random()

        getattr(solver, attr)(cval)


@pytest.mark.skipif(
    not amici.hdf5_enabled, reason="AMICI was compiled without HDF5"
)
def test_solver_hdf5_roundtrip(sbml_example_presimulation_module):
    """TestCase class for AMICI HDF5 I/O"""

    model = sbml_example_presimulation_module.get_model()
    solver = model.create_solver()
    _modify_solver_attrs(solver)

    hdf5file = "solverSettings.hdf5"

    amici.write_solver_settings_to_hdf5(solver, hdf5file, "ssettings")

    new_solver = model.create_solver()

    # check that we changed everything
    for attr in dir(solver):
        if not attr.startswith("set"):
            continue

        assert (
            getattr(solver, attr.replace("set", "get"))()
            != getattr(new_solver, attr.replace("set", "get"))()
        ), attr

    amici.read_solver_settings_from_hdf5(hdf5file, new_solver, "ssettings")

    # check that reading in settings worked
    for attr in dir(solver):
        if not attr.startswith("set"):
            continue

        assert getattr(solver, attr.replace("set", "get"))() == pytest.approx(
            getattr(new_solver, attr.replace("set", "get"))()
        ), attr

    os.remove(hdf5file)
