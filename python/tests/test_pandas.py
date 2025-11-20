"""Tests for ``amici.pandas``"""

import itertools

import numpy as np
import pytest
from amici.sim.sundials import (
    ExpData,
    get_data_observables_as_data_frame,
    get_edata_from_data_frame,
    run_simulation,
)
from amici.testing import skip_on_valgrind

# test parameters for test_pandas_import_export
combos = itertools.product([(10, 5), (5, 10), ()], repeat=3)
cases = [
    {
        "fixed_parameters": combo[0],
        "fixed_parameters_pre_equilibration": combo[1],
        "fixed_parameters_presimulation": combo[2],
    }
    for combo in combos
]


@skip_on_valgrind
@pytest.mark.parametrize("case", cases)
def test_pandas_import_export(sbml_example_presimulation_module, case):
    """TestCase class for testing csv import using pandas"""

    # setup
    model = sbml_example_presimulation_module.get_model()
    model.set_timepoints(np.linspace(0, 60, 61))
    solver = model.create_solver()
    rdata = run_simulation(model, solver)
    edata = [ExpData(rdata, 0.01, 0)]

    # test copy constructor
    _ = ExpData(edata[0])

    for fp in case:
        setattr(edata[0], fp, case[fp])

    df_edata = get_data_observables_as_data_frame(model, edata)
    edata_reconstructed = get_edata_from_data_frame(model, df_edata)

    for fp in [
        "fixed_parameters",
        "fixed_parameters_pre_equilibration",
        "fixed_parameters_presimulation",
    ]:
        if fp != "fixed_parameters" or case[fp] != ():
            assert getattr(edata[0], fp) == getattr(edata_reconstructed[0], fp)

            assert case[fp] == getattr(edata_reconstructed[0], fp)

        else:
            assert model.get_fixed_parameters() == getattr(
                edata_reconstructed[0], fp
            )

            assert model.get_fixed_parameters() == getattr(
                edata_reconstructed[0], fp
            )

        assert getattr(edata[0], fp) == case[fp]
