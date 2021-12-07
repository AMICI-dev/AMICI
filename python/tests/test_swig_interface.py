"""Test SWIG interface using python/examples/test_model_presimulation_pysb

Test getters, setters, etc.
"""

import copy
import numbers

import amici


def test_version_number(pysb_example_presimulation_module):
    model = pysb_example_presimulation_module.getModel()
    assert model.getAmiciVersion() == amici.__version__
    assert model.getAmiciCommit() == amici.__commit__


def test_copy_constructors(pysb_example_presimulation_module):
    model = pysb_example_presimulation_module.getModel()
    solver = model.getSolver()

    for obj in [model, solver]:
        for attr in dir(obj):
            if attr.startswith('__') \
                    or attr == 'this' \
                    or attr == 'thisown' \
                    or is_callable_but_not_getter(obj, attr):
                continue

            # objects will be initialized with default values so we
            # can check if the clone routine does the right thing by
            # modifying the value from default, cloning and checking
            # if the change was carried over to the clone
            val = get_val(obj, attr)

            try:
                modval = get_mod_val(val, attr)
            except ValueError:
                # happens for everything that is not bool or scalar
                continue

            try:
                set_val(obj, attr, modval)
            except AttributeError:
                # some attributes cannot be set
                continue

            obj_clone = obj.clone()

            assert get_val(obj, attr) == get_val(obj_clone, attr), \
                f"{obj} - {attr}"


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
    'AddSigmaResiduals': [
        False,
        True
    ],
    'AlwaysCheckFinite': [
        False,
        True,
    ],
    # Skipped due to model dependency in `'InitialStates'`.
    'FixedParameters': None,
    'InitialStates': [
        (10.0, 9.0, 1.0, 0.0, 0.0, 0.0),
        tuple([.1]*6),
    ],
    'InitialStateSensitivities': [
        tuple([1.0] + [0.0]*35),
        tuple([.1]*36),
    ],
    'MinimumSigmaResiduals': [
        50.0,
        60.0,
    ],
    ('nMaxEvent', 'setNMaxEvent'): [
        10,
        20,
    ],
    'Parameters': [
        (10.0, 0.1, 0.1, 0.1, 0.1, 0.1),
        tuple([1.0] * 6)
    ],
    # Skipped due to interdependency with `'InitialStateSensitivities'`.
    'ParameterList': None,
    # Skipped due to interdependency with `'InitialStateSensitivities'`.
    'ParameterScale': None,
    # Skipped due to interdependencies with
    # `'ReinitializeFixedParameterInitialStates'`.
    'ReinitializationStateIdxs': None,
    # Skipped due to interdependencies with `'ReinitializationStateIdxs'`.
    'ReinitializeFixedParameterInitialStates': None,
    # Skipped due to conservation laws in the test model
    # `pysb_example_presimulation_module.getModel()`.
    'StateIsNonNegative': None,
    'SteadyStateSensitivityMode': [
        0,
        1,
    ],
    ('t0', 'setT0'): [
        0.0,
        1.0,
    ],
    'Timepoints': [
        tuple(),
        (1.0, 2.0, 3.0),
    ],
}


def test_model_instance_settings(pysb_example_presimulation_module):
    model0 = pysb_example_presimulation_module.getModel()

    # Indexes of values in the `model_instance_settings0` dictionary.
    i_default = 0
    i_custom = 1

    i_getter = 0
    i_setter = 1

    # All settings are tested.
    assert set(model_instance_settings0) == set(amici.model_instance_settings)

    # Skip settings with interdependencies.
    model_instance_settings = \
        {k: v for k, v in model_instance_settings0.items() if v is not None}

    # All custom values are different to default values.
    assert all([
        default != custom
        for name, (default, custom) in model_instance_settings.items()
        if name != 'ReinitializeFixedParameterInitialStates'
    ])

    # All default values are as expected.
    for name, (default, custom) in model_instance_settings.items():
        getter = name[i_getter] if isinstance(name, tuple) else f'get{name}'
        setter = name[i_setter] if isinstance(name, tuple) else f'set{name}'
        # Default values are as expected.
        assert getattr(model0, getter)() == default
        # Custom value is set correctly.
        getattr(model0, setter)(custom)
        assert getattr(model0, getter)() == custom

    # The grouped getter method works.
    custom_settings = amici.get_model_settings(model0)
    for name in model_instance_settings:
        assert custom_settings[name] == model_instance_settings[name][i_custom]

    # Create a new model for comparison.
    model = pysb_example_presimulation_module.getModel()

    # The new model has the default settings.
    model_default_settings = amici.get_model_settings(model)
    for name in model_instance_settings:
        assert model_default_settings[name] == \
            model_instance_settings[name][i_default]

    # The grouped setter method works.
    custom_settings_not_none = {
        name: value
        for name, value in custom_settings.items()
        if model_instance_settings0[name] is not None
    }
    amici.set_model_settings(model, custom_settings_not_none)
    assert all([
        value == custom_settings_not_none[name]
        for name, value in amici.get_model_settings(model).items()
        if name in custom_settings_not_none
    ])


def test_interdependent_settings(pysb_example_presimulation_module):
    """Test settings that were not tested in `test_model_instance_settings`.

    `StateIsNonNegative` is still skipped, due to conservation laws in the
    test model.
    """
    model = pysb_example_presimulation_module.getModel()

    original_settings = {
        'FixedParameters': (9.0, 1.0),
        'ParameterList': (0, 1, 2, 3, 4, 5),
        'ParameterScale': [0, 0, 0, 0, 0, 0],
        'ReinitializationStateIdxs': tuple(),
        'ReinitializeFixedParameterInitialStates': False,
        'StateIsNonNegative': (False, False, False),
    }

    expected_settings = {
        'FixedParameters': (8.0, 2.0),
        'ParameterList': (0, 1, 2, 3, 4),
        'ParameterScale': [1, 0, 0, 0, 0, 0],
        'ReinitializationStateIdxs': (0,),
        'ReinitializeFixedParameterInitialStates': True,
        # Skipped due to conservation laws in the test model.
        # 'StateIsNonNegative': None,
    }

    ignored_settings = [
        # Ignored due to conservation laws in the test model.
        'StateIsNonNegative',
    ]

    # Some values need to be transformed to be tested in Python
    # (e.g. SWIG objects). Default transformer is no transformation
    # (the identity function).
    getter_transformers = {
        setting: (lambda x: x)
        for setting in expected_settings
    }
    getter_transformers.update({
        # Convert from SWIG object.
        'ParameterScale': lambda x: list(x)
    })

    default_settings = {
        setting: setting_value
        for setting, setting_value in amici.get_model_settings(model).items()
        if setting not in ignored_settings
    }

    for setting, expected_value in expected_settings.items():
        input_settings = {setting: copy.deepcopy(expected_value)}

        amici.set_model_settings(model, input_settings)
        output_settings = amici.get_model_settings(model)
        test_value = getter_transformers[setting](
            output_settings[setting]
        )
        # The setter works.
        assert test_value == expected_value

        input_settings = {setting: output_settings[setting]}
        amici.set_model_settings(model, input_settings)
        output_settings = amici.get_model_settings(model)
        test_value = getter_transformers[setting](
            output_settings[setting]
        )
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
    model = pysb_example_presimulation_module.getModel()

    not_handled = [
        'get',
        'getAmiciCommit',
        'getAmiciVersion',
        'getExpressionIds',
        'getExpressionNames',
        'getFixedParameterById',
        'getFixedParameterByName',
        'getFixedParameterIds',
        'getFixedParameterNames',
        'getName',
        'getObservableIds',
        'getObservableNames',
        'getObservableScaling',
        'getParameterById',
        'getParameterByName',
        'getParameterIds',
        'getParameterNames',
        'getSolver',
        'getStateIds',
        'getStateNames',
        'getTimepoint',
        'getUnscaledParameters',
        'setAllStatesNonNegative',
        'setFixedParameterById',
        'setFixedParameterByName',
        'setFixedParametersByIdRegex',
        'setFixedParametersByNameRegex',
        'setParameterById',
        'setParameterByName',
        'setParametersByIdRegex',
        'setParametersByNameRegex',
        'setUnscaledInitialStateSensitivities',
    ]

    handled = [
        name
        for names in amici.model_instance_settings
        for name in (
            names
            if isinstance(names, tuple) else
            (f'get{names}', f'set{names}')
        )
    ]

    for attribute in dir(model):
        if attribute[:3] in ['get', 'set'] and attribute not in not_handled:
            assert attribute in handled, attribute


def is_callable_but_not_getter(obj, attr):
    if not callable(getattr(obj, attr)):
        return False

    if attr.startswith('get'):
        return \
            'set' + attr[3:] not in dir(obj) \
            or attr.endswith('ById') \
            or attr.endswith('ByName')
    else:
        return True


def get_val(obj, attr):
    if callable(getattr(obj, attr)):
        return getattr(obj, attr)()
    else:
        return getattr(obj, attr)


def get_mod_val(val, attr):
    if attr == 'getReturnDataReportingMode':
        return amici.RDataReporting.likelihood
    elif attr == 'getParameterList':
        return tuple(get_mod_val(val[0], '') for _ in val)
    elif attr == 'getStateIsNonNegative':
        raise ValueError('Cannot modify value')
    elif isinstance(val, bool):
        return not val
    elif isinstance(val, numbers.Number):
        return val + 1
    elif isinstance(val, tuple):
        return tuple(get_mod_val(v, attr) for v in val)

    raise ValueError('Cannot modify value')


def set_val(obj, attr, val):
    if callable(getattr(obj, attr)):
        getattr(obj, 'set' + attr[3:])(
            val
        )
    else:
        setattr(obj, attr, val)
