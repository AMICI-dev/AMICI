"""Test SWIG interface using python/examples/test_model_presimulation_pysb

Test getters, setters, etc.
"""

import copy
import numbers
import pickle
from math import nan

import amici
import numpy as np
import pytest
import xarray
from amici.sim.sundials import (
    AMICI_FIRST_RHSFUNC_ERR,
    CVodeSolver,
    ExpData,
    ExpDataPtr,
    ExpDataView,
    IDASolver,
    ModelPtr,
    ParameterScaling,
    RDataReporting,
    ReturnDataView,
    SensitivityMethod,
    SensitivityOrder,
    SolverPtr,
    SteadyStateComputationMode,
    SteadyStateSensitivityMode,
    get_model_settings,
    parameter_scaling_from_int_vector,
    run_simulation,
    run_simulations,
    set_model_settings,
)
from amici.testing import skip_on_valgrind


def test_version_number(pysb_example_presimulation_module):
    model = pysb_example_presimulation_module.get_model()
    assert model.get_amici_version() == amici.__version__
    assert model.get_amici_commit() == amici.__commit__


def test_copy_constructors(pysb_example_presimulation_module):
    model = pysb_example_presimulation_module.get_model()
    solver = model.create_solver()

    for obj in [model, solver]:
        for attr in dir(obj):
            if (
                attr.startswith("__")
                or attr == "this"
                or attr == "thisown"
                or is_callable_but_not_getter(obj, attr)
            ):
                continue

            # objects will be initialized with default values so we
            # can check if the clone routine does the right thing by
            # modifying the value from default, cloning and checking
            # if the change was carried over to the clone
            val = get_val(obj, attr)

            try:
                modval = get_mod_val(val, attr, obj)
            except ValueError:
                # happens for everything that is not bool or scalar
                continue

            try:
                set_val(obj, attr, modval)
            except AttributeError:
                # some attributes cannot be set
                continue

            obj_clone = obj.clone()

            assert get_val(obj, attr) == get_val(obj_clone, attr), (
                f"{obj} - {attr}"
            )


# `None` values are skipped in `test_model_instance_settings`.
# Keys are suffixes of `get[...]` and `set[...]` `amici.Model` methods.
# If either the getter or setter is not named with this pattern, then the key
# is a tuple where the first and second elements are the getter and setter
# methods, respectively.
# Values are lists where the first element is the default value of the test
# model, and the second value is some custom value.
# Default values are based on `pysb_example_presimulation_module`.
model_instance_settings0 = {
    # setting name: [default value, custom value]
    "add_sigma_residuals": [False, True],
    "always_check_finite": [False, True],
    # Skipped due to model dependency in `'InitialStates'`.
    "fixed_parameters": None,
    "initial_state": [
        (10.0, 9.0, 1.0, 0.0, 0.0, 0.0),
        tuple([0.1] * 6),
    ],
    (
        "get_initial_state_sensitivities",
        "set_unscaled_initial_state_sensitivities",
    ): [
        tuple([1.0] + [0.0] * 35),
        tuple([0.1] * 36),
    ],
    "steadystate_mask": [
        (),
        tuple([0] * 3),
    ],
    "minimum_sigma_residuals": [
        50.0,
        60.0,
    ],
    ("n_max_event", "set_n_max_event"): [
        10,
        20,
    ],
    "free_parameters": [(10.0, 0.1, 0.1, 0.1, 0.1, 0.1), tuple([1.0] * 6)],
    # Skipped due to interdependency with `'InitialStateSensitivities'`.
    "parameter_list": None,
    # Skipped due to interdependency with `'InitialStateSensitivities'`.
    "parameter_scale": None,
    # Skipped due to interdependencies with
    # `'ReinitializeFixedParameterInitialStates'`.
    "reinitialization_state_idxs": None,
    # Skipped due to interdependencies with `'ReinitializationStateIdxs'`.
    "reinitialize_fixed_parameter_initial_states": None,
    # Skipped due to conservation laws in the test model
    # `pysb_example_presimulation_module.getModel()`.
    "state_is_non_negative": None,
    "steady_state_computation_mode": [
        SteadyStateComputationMode.integrationOnly,
        SteadyStateComputationMode.integrateIfNewtonFails,
    ],
    "steady_state_sensitivity_mode": [
        SteadyStateSensitivityMode.integrationOnly,
        SteadyStateSensitivityMode.integrateIfNewtonFails,
    ],
    ("t0", "set_t0"): [
        0.0,
        1.0,
    ],
    ("t0_preeq", "set_t0_preeq"): [
        nan,
        -10.0,
    ],
    "timepoints": [
        tuple(),
        (1.0, 2.0, 3.0),
    ],
}


def same_or_nan(a, b):
    """Check if two values are the same or both NaN."""
    return a == b or (isinstance(a, float) and np.isnan(a) and np.isnan(b))


def test_model_instance_settings(pysb_example_presimulation_module):
    model0 = pysb_example_presimulation_module.get_model()

    # Indexes of values in the `model_instance_settings0` dictionary.
    i_default = 0
    i_custom = 1

    i_getter = 0
    i_setter = 1

    # the default setting for AlwaysCheckFinite depends on whether the amici
    # extension has been built in debug mode
    model_instance_settings0["always_check_finite"] = [
        model0.get_always_check_finite(),
        not model0.get_always_check_finite(),
    ]

    # All settings are tested.
    assert set(model_instance_settings0) == set(
        amici.sim.sundials._swig_wrappers.model_instance_settings
    )

    # Skip settings with interdependencies.
    model_instance_settings = {
        k: v for k, v in model_instance_settings0.items() if v is not None
    }

    # All custom values are different to default values.
    assert all(
        default != custom
        for name, (default, custom) in model_instance_settings.items()
        if name != "reinitialize_fixed_parameter_initial_states"
    )

    # All default values are as expected.
    for name, (default, custom) in model_instance_settings.items():
        getter = name[i_getter] if isinstance(name, tuple) else f"get_{name}"
        setter = name[i_setter] if isinstance(name, tuple) else f"set_{name}"
        # Default values are as expected.
        assert same_or_nan(getattr(model0, getter)(), default), name
        # Custom value is set correctly.
        getattr(model0, setter)(custom)
        assert getattr(model0, getter)() == custom

    # The grouped getter method works.
    custom_settings = get_model_settings(model0)
    for name in model_instance_settings:
        assert custom_settings[name] == model_instance_settings[name][i_custom]

    # Create a new model for comparison.
    model = pysb_example_presimulation_module.get_model()

    # The new model has the default settings.
    model_default_settings = get_model_settings(model)
    for name in model_instance_settings:
        if (
            name == "initial_state" and not model.has_custom_initial_state()
        ) or (
            name
            == (
                "get_initial_state_sensitivities",
                "set_unscaled_initial_state_sensitivities",
            )
            and not model.has_custom_initial_state_sensitivities()
        ):
            # Here the expected value differs from what the getter would return
            assert model_default_settings[name] == []
        else:
            assert same_or_nan(
                model_default_settings[name],
                model_instance_settings[name][i_default],
            ), name

    # The grouped setter method works.
    custom_settings_not_none = {
        name: value
        for name, value in custom_settings.items()
        if model_instance_settings0[name] is not None
    }
    set_model_settings(model, custom_settings_not_none)
    assert all(
        same_or_nan(value, custom_settings_not_none[name])
        for name, value in get_model_settings(model).items()
        if name in custom_settings_not_none
    )


def test_interdependent_settings(pysb_example_presimulation_module):
    """Test settings that were not tested in `test_model_instance_settings`.

    `StateIsNonNegative` is still skipped, due to conservation laws in the
    test model.
    """
    model = pysb_example_presimulation_module.get_model()

    original_settings = {
        "fixed_parameters": (9.0, 1.0),
        "parameter_list": (0, 1, 2, 3, 4, 5),
        "parameter_scale": [0, 0, 0, 0, 0, 0],
        "reinitialization_state_idxs": tuple(),
        "reinitialize_fixed_parameter_initial_states": False,
        "state_is_non_negative": (False, False, False),
    }

    expected_settings = {
        "fixed_parameters": (8.0, 2.0),
        "parameter_list": (0, 1, 2, 3, 4),
        "parameter_scale": [1, 0, 0, 0, 0, 0],
        "reinitialization_state_idxs": (0,),
        "reinitialize_fixed_parameter_initial_states": True,
        # Skipped due to conservation laws in the test model.
        # 'state_is_non_negative': None,
    }

    # Some values need to be transformed to be tested in Python
    # (e.g. SWIG objects). Default transformer is no transformation
    # (the identity function).
    getter_transformers = {
        setting: (lambda x: x) for setting in original_settings
    }
    getter_transformers.update(
        {
            # Convert from SWIG object.
            "parameter_scale": lambda x: list(x)
        }
    )

    default_settings = get_model_settings(model)
    for original_setting, original_setting_value in original_settings.items():
        test_value = getter_transformers[original_setting](
            default_settings[original_setting],
        )
        # The original model is as expected (to ensure the test is still
        # valid).
        assert test_value == original_setting_value

    for setting, expected_value in expected_settings.items():
        input_settings = {setting: copy.deepcopy(expected_value)}

        set_model_settings(model, input_settings)
        output_settings = get_model_settings(model)
        test_value = getter_transformers[setting](output_settings[setting])
        # The setter works.
        assert test_value == expected_value

        input_settings = {setting: output_settings[setting]}
        set_model_settings(model, input_settings)
        output_settings = get_model_settings(model)
        test_value = getter_transformers[setting](output_settings[setting])
        # (round-trip) The output of the getter can be used as input to the
        # setter, and does not change the value.
        assert test_value == expected_value


def test_unhandled_settings(pysb_example_presimulation_module):
    """Detect possible getters and setters that are not yet handled.

    Currently, only model attributes that begin with `'get'` and `'set'` are
    tested.

    The following getters/setters are untested here, due to unusual naming
    (not prefixed by `'get'` or `'set'`).
    - nMaxEvent
    - t0
    """
    model = pysb_example_presimulation_module.get_model()

    not_handled = [
        "get",
        "get_amici_commit",
        "get_amici_version",
        "get_explicit_roots",
        "get_expression_ids",
        "get_expression_names",
        "get_fixed_parameter_by_id",
        "get_fixed_parameter_by_name",
        "get_fixed_parameter_ids",
        "get_fixed_parameter_names",
        "get_id_list",
        "get_name",
        "get_observable_ids",
        "get_observable_names",
        "get_observable_scaling",
        "get_free_parameter_by_id",
        "get_free_parameter_by_name",
        "get_free_parameter_ids",
        "get_free_parameter_names",
        "get_second_order_mode",
        "get_solver",
        "get_state_ids",
        "get_state_names",
        "get_state_ids_solver",
        "get_state_names_solver",
        "get_timepoint",
        "get_unscaled_parameters",
        "set_all_states_non_negative",
        "set_fixed_parameter_by_id",
        "set_fixed_parameter_by_name",
        "set_fixed_parameters_by_id_regex",
        "set_fixed_parameters_by_name_regex",
        "set_free_parameter_by_id",
        "set_free_parameter_by_name",
        "set_free_parameters_by_id_regex",
        "set_free_parameters_by_name_regex",
        "set_initial_state_sensitivities",
        "get_trigger_timepoints",
        "get_any_state_nonnegative",
    ]
    from amici.sim.sundials._swig_wrappers import model_instance_settings

    handled = [
        name
        for names in model_instance_settings
        for name in (
            names
            if isinstance(names, tuple)
            else (f"get_{names}", f"set_{names}")
        )
    ]

    for attribute in dir(model):
        if attribute[:3] in ["get", "set"] and attribute not in not_handled:
            assert attribute in handled, attribute


def is_callable_but_not_getter(obj, attr):
    if not callable(getattr(obj, attr)):
        return False

    if attr.startswith("get"):
        return (
            "set" + attr[3:] not in dir(obj)
            or attr.endswith("by_id")
            or attr.endswith("by_name")
        )
    else:
        return True


def get_val(obj, attr):
    if callable(getattr(obj, attr)):
        return getattr(obj, attr)()
    else:
        return getattr(obj, attr)


def get_mod_val(val, attr, obj):
    if attr == "get_return_data_reporting_mode":
        return RDataReporting.likelihood
    elif attr == "get_parameter_list":
        return tuple(get_mod_val(val[0], "", obj) for _ in val)
    elif attr == "get_state_is_non_negative":
        raise ValueError("Cannot modify value")
    elif attr == "get_steadystate_mask":
        return [0 for _ in range(obj.nx_solver)]
    elif isinstance(val, bool):
        return not val
    elif isinstance(val, numbers.Number):
        return val + 1
    elif isinstance(val, tuple):
        return tuple(get_mod_val(v, attr, obj) for v in val)

    raise ValueError("Cannot modify value")


def set_val(obj, attr, val):
    if callable(getattr(obj, attr)):
        getattr(obj, "set" + attr[3:])(val)
    else:
        setattr(obj, attr, val)


def test_model_instance_settings_custom_x0(pysb_example_presimulation_module):
    """Check that settings are applied in the correct order, and only if
    required"""

    def assert_same(a: dict, b: dict):
        """Assert that two model settings dictionaries are the same."""
        assert set(a.keys()) == set(b.keys()), a.keys() ^ b.keys()
        for key in a:
            assert same_or_nan(a[key], b[key]), f"{key}: {a[key]} != {b[key]}"

    model = pysb_example_presimulation_module.get_model()

    # ensure no-custom-(s)x0 is restored
    assert not model.has_custom_initial_state()
    assert not model.has_custom_initial_state_sensitivities()
    settings = get_model_settings(model)
    model.set_initial_state(model.get_initial_state())
    model.set_unscaled_initial_state_sensitivities(
        model.get_initial_state_sensitivities()
    )
    set_model_settings(model, settings)
    assert not model.has_custom_initial_state()
    assert not model.has_custom_initial_state_sensitivities()
    # ensure everything was set correctly, and there wasn't any problem
    #  due to, e.g., interactions of different setters
    assert_same(settings, get_model_settings(model))

    # ensure custom (s)x0 is restored
    model.set_initial_state(model.get_initial_state())
    model.set_parameter_scale(ParameterScaling.log10)
    sx0 = model.get_initial_state_sensitivities()
    model.set_unscaled_initial_state_sensitivities(sx0)
    assert model.has_custom_initial_state()
    assert model.has_custom_initial_state_sensitivities()
    settings = get_model_settings(model)
    model2 = pysb_example_presimulation_module.get_model()
    set_model_settings(model2, settings)
    assert model2.has_custom_initial_state()
    assert model2.has_custom_initial_state_sensitivities()
    assert model2.get_initial_state_sensitivities() == sx0
    assert_same(settings, get_model_settings(model2))


def test_solver_repr():
    for solver in (CVodeSolver(), IDASolver()):
        solver_ptr = SolverPtr(solver.this)
        for s in (solver, solver_ptr):
            assert "maxsteps" in str(s)
            assert "maxsteps" in repr(s)


def test_edata_repr():
    ny = 1
    nz = 2
    ne = 3
    nt = 4
    edata = ExpData(ny, nz, ne, range(nt))
    edata_ptr = ExpDataPtr(edata.this)
    expected_strs = (
        f"{nt}x{ny} time-resolved datapoints",
        f"{ne}x{nz} event-resolved datapoints",
        f"(0/{ny * nt} measurements",
        f"(0/{nz * ne} measurements",
    )
    for e in [edata, edata_ptr]:
        for expected_str in expected_strs:
            assert expected_str in str(e)
            assert expected_str in repr(e)


def test_edata_equality_operator():
    e1 = ExpData(1, 2, 3, [3])
    e2 = ExpData(1, 2, 3, [3])
    assert e1 == e2
    # check that comparison with other types works
    # this is not implemented by swig by default
    assert e1 != 1


def test_expdata_and_expdataview_are_deepcopyable():
    edata1 = ExpData(3, 2, 3, range(4))
    edata1.set_measurements(np.zeros((3, 4)).flatten())

    # ExpData
    edata2 = copy.deepcopy(edata1)
    assert edata1 == edata2
    assert edata1.this != edata2.this
    edata2.set_timepoints([0])
    assert edata1 != edata2

    # ExpDataView
    ev1 = ExpDataView(edata1)
    ev2 = copy.deepcopy(ev1)
    assert ev2._swigptr.this != ev1._swigptr.this
    assert ev1 == ev2


def test_solvers_are_deepcopyable():
    for solver_type in (CVodeSolver, IDASolver):
        for solver1 in (solver_type(), SolverPtr(solver_type())):
            solver2 = copy.deepcopy(solver1)
            assert solver1.this != solver2.this
            assert (
                solver1.get_relative_tolerance()
                == solver2.get_relative_tolerance()
            )
            solver2.set_relative_tolerance(
                100 * solver2.get_relative_tolerance()
            )
            assert (
                solver1.get_relative_tolerance()
                != solver2.get_relative_tolerance()
            )


def test_model_is_deepcopyable(pysb_example_presimulation_module):
    model_module = pysb_example_presimulation_module
    for model1 in (
        model_module.get_model(),
        ModelPtr(model_module.get_model()),
    ):
        model2 = copy.deepcopy(model1)
        assert model1.this != model2.this
        assert model1.t0() == model2.t0()
        model2.set_t0(100 + model2.t0())
        assert model1.t0() != model2.t0()


def test_rdataview(sbml_example_presimulation_module):
    """Test some SwigPtrView functionality via ReturnDataView."""
    model_module = sbml_example_presimulation_module
    model = model_module.get_model()
    model.set_timepoints([1, 2, 3])
    rdata = run_simulation(model, model.create_solver())
    assert isinstance(rdata, ReturnDataView)

    # check that non-array attributes are looked up in the wrapped object
    assert rdata.ptr.ny == rdata.ny

    # fields are accessible via dot notation and [] operator,
    #  __contains__ and __getattr__ are implemented correctly
    with pytest.raises(AttributeError):
        _ = rdata.nonexisting_attribute

    with pytest.raises(KeyError):
        _ = rdata["nonexisting_attribute"]

    assert not hasattr(rdata, "nonexisting_attribute")
    assert "x" in rdata
    assert (rdata.x == rdata["x"]).all()

    # field names are included by dir()
    assert "x" in dir(rdata)

    # Test xarray conversion
    xr_x = rdata.xr.x
    assert isinstance(xr_x, xarray.DataArray)
    assert (rdata.x == xr_x).all()
    assert xr_x.dims == ("time", "state")
    assert (xr_x.coords["time"].data == rdata.ts).all()
    assert (xr_x.coords["state"].data == model.get_state_ids()).all()

    # test that generating the xarrays does not fail, without checking
    #  their content any further
    for attr in dir(rdata.xr):
        if not attr.startswith("_"):
            try:
                getattr(rdata.xr, attr)
            except TypeError as e:
                print(str(e))


def test_rdata_ids(sbml_example_presimulation_module):
    """Test that rdata IDs are correctly set."""
    model_module = sbml_example_presimulation_module
    model = model_module.get_model()

    model.set_timepoints([0, 1, 2])
    rdata = model.simulate()

    assert isinstance(rdata.free_parameter_ids, tuple)
    assert rdata.free_parameter_ids == model.get_free_parameter_ids()
    assert rdata.fixed_parameter_ids == model.get_fixed_parameter_ids()
    assert rdata.state_ids == model.get_state_ids()
    assert rdata.state_ids_solver == model.get_state_ids_solver()
    assert rdata.observable_ids == model.get_observable_ids()
    assert rdata.expression_ids == model.get_expression_ids()


def test_python_exceptions(sbml_example_presimulation_module):
    """Test that C++ exceptions are correctly caught and re-raised in Python."""
    from amici.sim.sundials import run_simulation

    # amici-base extension throws and its swig-wrapper catches
    solver = CVodeSolver()
    with pytest.raises(
        RuntimeError, match="maxsteps must be a positive number"
    ):
        solver.set_max_steps(-1)

    # model extension throws and its swig-wrapper catches
    model = sbml_example_presimulation_module.get_model()
    with pytest.raises(RuntimeError, match="Steadystate mask has wrong size"):
        model.set_steadystate_mask([1] * model.nx_solver * 2)

    # amici-base extension throws and its swig-wrapper catches
    edata = ExpData(1, 1, 1, [1])
    # too short sx0
    edata.sx0 = (1, 2)
    with pytest.raises(
        RuntimeError,
        match=r"Number of initial conditions sensitivities \(36\) "
        r"in model does not match ExpData \(2\).",
    ):
        run_simulation(model, solver, edata)

    run_simulations(
        model, solver, [edata, edata], failfast=True, num_threads=1
    )

    # model throws, base catches, swig-exception handling is not involved
    model.set_free_parameters([nan] * model.np())
    model.set_timepoints([1])
    rdata = run_simulation(model, solver)
    assert rdata.status == AMICI_FIRST_RHSFUNC_ERR

    edata = ExpData(1, 1, 1, [1])
    rdatas = run_simulations(
        model, solver, [edata, edata], failfast=True, num_threads=1
    )
    assert rdatas[0].status == AMICI_FIRST_RHSFUNC_ERR

    # model throws, base catches, swig-exception handling is involved
    with pytest.raises(
        RuntimeError, match="AMICI failed to integrate the forward problem"
    ):
        # rethrow=True
        amici._installation._amici.run_simulation(
            solver, None, model.get(), True
        )


def test_reporting_mode_obs_llh(sbml_example_presimulation_module):
    model_module = sbml_example_presimulation_module
    model = model_module.get_model()
    solver = model.create_solver()

    solver.set_return_data_reporting_mode(
        RDataReporting.observables_likelihood
    )
    solver.set_sensitivity_order(SensitivityOrder.first)

    for sens_method in (
        SensitivityMethod.none,
        SensitivityMethod.forward,
        SensitivityMethod.adjoint,
    ):
        solver.set_sensitivity_method(sens_method)
        rdata = run_simulation(model, solver, ExpData(1, 1, 1, [1]))
        assert rdata.rdata_reporting == RDataReporting.observables_likelihood

        assert rdata.y.size > 0
        assert rdata.sigmay.size > 0
        assert rdata.J is None

        match solver.get_sensitivity_method():
            case SensitivityMethod.none:
                assert rdata.sllh is None
            case SensitivityMethod.forward:
                assert rdata.sy.size > 0
                assert rdata.ssigmay.size > 0
                assert rdata.sllh.size > 0
                assert not np.isnan(rdata.sllh).any()
            case SensitivityMethod.adjoint:
                assert rdata.sy is None
                assert rdata.ssigmay is None
                assert rdata.sllh.size > 0
                assert not np.isnan(rdata.sllh).any()


@skip_on_valgrind
def test_pickle_model(sbml_example_presimulation_module):
    model_module = sbml_example_presimulation_module
    model = model_module.get_model()

    assert (
        model.get_steady_state_sensitivity_mode()
        == SteadyStateSensitivityMode.integrationOnly
    )
    model.set_steady_state_sensitivity_mode(
        SteadyStateSensitivityMode.newtonOnly
    )

    model_pickled = pickle.loads(pickle.dumps(model))
    # ensure it's re-picklable
    model_pickled = pickle.loads(pickle.dumps(model_pickled))
    assert (
        model_pickled.get_steady_state_sensitivity_mode()
        == SteadyStateSensitivityMode.newtonOnly
    )

    model_pickled.set_steady_state_sensitivity_mode(
        SteadyStateSensitivityMode.integrateIfNewtonFails
    )
    assert (
        model.get_steady_state_sensitivity_mode()
        != model_pickled.get_steady_state_sensitivity_mode()
    )

    # ensure we can pickle after clone()
    model_clone = model.clone()
    pickle.loads(pickle.dumps(model_clone))


def test_pickle_edata():
    ny = 2
    nz = 3
    ne = 4
    nt = 5
    edata = ExpData(ny, nz, ne, range(nt))
    edata.set_measurements(list(np.arange(ny * nt, dtype=float)))
    edata.pscale = parameter_scaling_from_int_vector(
        [ParameterScaling.log10] * 5
    )

    edata_pickled = pickle.loads(pickle.dumps(edata))
    assert edata == edata_pickled


@pytest.mark.skipif(
    not amici.sim.sundials.hdf5_enabled,
    reason="AMICI build without HDF5 support",
)
def test_pickle_solver():
    for solver in (
        CVodeSolver(),
        IDASolver(),
        SolverPtr(CVodeSolver()),
        SolverPtr(IDASolver()),
    ):
        solver.set_max_steps(1234)
        solver.set_sensitivity_order(SensitivityOrder.first)
        solver_pickled = pickle.loads(pickle.dumps(solver))
        assert type(solver) is type(solver_pickled)
        assert solver.get_max_steps() == solver_pickled.get_max_steps()
        assert (
            solver.get_sensitivity_order()
            == solver_pickled.get_sensitivity_order()
        )
