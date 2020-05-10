"""AMICI HDF5 I/O tests"""

import os
import random

import amici
import pytest


def _modify_solver_attrs(solver):
    # change to non-default values
    for attr in dir(solver):
        if not attr.startswith('set'):
            continue

        val = getattr(solver, attr.replace('set', 'get'))()

        if isinstance(val, bool):
            cval = not val
        elif attr == 'setStabilityLimitFlag':
            cval = 0
        elif attr == 'setReturnDataReportingMode':
            cval = amici.RDataReporting.likelihood
        elif isinstance(val, int):
            cval = val + 1
        else:
            cval = val + random.random()

        getattr(solver, attr)(cval)


def test_solver_hdf5_roundtrip(sbml_example_presimulation_module):
    """TestCase class for AMICI HDF5 I/O"""

    model = sbml_example_presimulation_module.getModel()
    solver = model.getSolver()
    _modify_solver_attrs(solver)

    hdf5file = 'solverSettings.hdf5'

    amici.writeSolverSettingsToHDF5(solver, hdf5file, 'ssettings')

    new_solver = model.getSolver()

    # check that we changed everything
    for attr in dir(solver):
        if not attr.startswith('set'):
            continue

        assert getattr(solver, attr.replace('set', 'get'))() \
            != getattr(new_solver, attr.replace('set', 'get'))(), attr

    amici.readSolverSettingsFromHDF5(hdf5file, new_solver, 'ssettings')

    # check that reading in settings worked
    for attr in dir(solver):
        if not attr.startswith('set'):
            continue

        assert getattr(solver, attr.replace('set', 'get'))() \
            == pytest.approx(
                getattr(new_solver, attr.replace('set', 'get'))()), attr

    os.remove(hdf5file)
