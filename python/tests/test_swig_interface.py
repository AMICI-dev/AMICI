"""Test SWIG interface using python/examples/test_model_presimulation_pysb

Test getters, setters, etc.
"""

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
