import os
import tempfile
import math
import numpy as np
import sympy as sp
import pandas as pd

# from contextlib import contextmanager
from scipy.integrate import quad

import libsbml
import petab
import amici

from amici.petab_import import import_petab_problem
from amici.petab_objective import (
    simulate_petab,
    create_parameterized_edatas,
    LLH,
    SLLH,
    RDATAS,
)
from amici.gradient_check import check_derivatives, check_close, check_results
from amici.sbml_utils import (
    amici_time_symbol,
    createSbmlModel,
    addCompartment,
    addParameter,
    addSpecies,
    addRateRule,
    addInflow,
    setSbmlMath,
)
from amici.splines import AbstractSpline, CubicHermiteSpline, UniformGrid


# @contextmanager
# def cd(path):
#     cwd = os.getcwd()
#     assert os.path.isabs(cwd)
#     os.chdir(path)
#     yield os.path.abspath(path)
#     os.chdir(pwd)


def assert_fun(x):
    assert x


def evaluate_spline(spline, params, tt, **kwargs):
    return np.asarray([spline.evaluate(t).subs(params) for t in tt], **kwargs)


def integrate_spline(spline, params, tt, initial_value=0, **kwargs):
    ispline = [initial_value + spline.integrate(0, t) for t in tt]
    if params is not None:
        ispline = [x.subs(params) for x in ispline]
    return np.asarray(ispline, **kwargs)


def create_condition_table():
    condition_df = pd.DataFrame({'conditionId' : ['condition1']})
    condition_df.set_index(['conditionId'], inplace=True)
    return condition_df


def create_parameter_table(**kwargs):
    if isinstance(kwargs['parameterId'], str):
        kwargs['parameterId'] = [kwargs['parameterId']]
    if 'parameterScale' not in kwargs.keys():
        kwargs['parameterScale'] = 'lin'
    if 'estimate' not in kwargs.keys():
        kwargs['estimate'] = 1
    parameter_df = pd.DataFrame(kwargs)
    parameter_df.set_index(['parameterId'], inplace=True)
    return parameter_df


def create_observable_table(**kwargs):
    if isinstance(kwargs['observableId'], str):
        kwargs['observableId'] = [kwargs['observableId']]
    if 'observableTransformation' not in kwargs.keys():
        kwargs['observableTransformation'] = 'lin'
    if 'noiseDistribution' not in kwargs.keys():
        kwargs['noiseDistribution'] = 'normal'
    observable_df = pd.DataFrame(kwargs)
    observable_df.set_index(['observableId'], inplace=True)
    return observable_df


def create_measurement_table(**kwargs):
    if isinstance(kwargs['observableId'], str):
        kwargs['observableId'] = [kwargs['observableId']]
    if 'simulationConditionId' not in kwargs.keys():
        kwargs['simulationConditionId'] = 'condition1'
    return pd.DataFrame(kwargs)


def species(i):
    return f'z{i}'


def observable(i):
    return f'{species(i)}_obs'


def species_to_index(name):
    assert name[0] == 'z'
    return int(name[1:])


def create_petab_problem(splines,
                         params_true,
                         initial_values,
                         use_reactions=False,
                         measure_upsample=3,
                         sigma=1.0,
                         Textrapolate=0.25,
                         folder=None
    ):

    modelname = 'test_splines'

    for spline in splines:
        assert spline.x == amici_time_symbol

    # Create SBML document
    doc, model = createSbmlModel(modelname)
    addCompartment(model, 'compartment')
    for i in range(len(splines)):
        splines[i].addToSbmlModel(model)
        addSpecies(model, species(i), initial_amount=initial_values[i])
        if use_reactions:
            addInflow(model, species(i), splines[i].sbmlId)
        else:
            addRateRule(model, species(i), splines[i].sbmlId)
    for (parId, value) in params_true.items():
        addParameter(model, parId, value=value, constant=True)
    for spline in splines:
        addParameter(model, spline.sbmlId, constant=False)

    # Compute simulation time
    T = 0
    for spline in splines:
        if spline.extrapolate[0] is None:
            assert spline.xx[0] <= 0
        if spline.extrapolate[1] is not None:
            DT = (spline.xx[-1] - spline.xx[0]) * Textrapolate
        else:
            DT = 0
        T = max(T, spline.xx[-1] + DT)

    # Compute measurements
    dt = min(np.diff(spline.xx).min() for spline in splines)
    dt /= measure_upsample
    n_obs = math.ceil(T / dt) + 1
    tt_obs = np.linspace(0, float(T), n_obs)
    zz_true = [
        integrate_spline(spline, params_true, tt_obs, iv, dtype=float)
        for (spline, iv) in zip(splines, initial_values)
    ]
    zz_obs = [zz + sigma * np.random.randn(len(zz)) for zz in zz_true]

    # Create PEtab tables
    condition_df = create_condition_table()
    _params = list(params_true.items())
    parameter_df = create_parameter_table(
        parameterId = [p.name for (p, v) in _params],
        lowerBound = min(v for (p, v) in _params),
        upperBound = max(v for (p, v) in _params),
        nominalValue = [v for (p, v) in _params],
        estimate = 1
    )
    observable_df = create_observable_table(
        observableId = [observable(i) for i in range(len(splines))],
        observableFormula = [species(i) for i in range(len(splines))],
        noiseFormula = sigma if sigma > 0 else 1.0
    )
    measurement_df = create_measurement_table(
        observableId = np.concatenate([len(tt_obs) * [observable(i)] for i in range(len(splines))]),
        time = len(splines) * list(tt_obs),
        measurement = np.concatenate(zz_obs)
    )

    # Create and validate PEtab problem
    problem = petab.Problem(
        sbml_document = doc,
        sbml_model = model,
        condition_df = condition_df,
        measurement_df = measurement_df,
        parameter_df = parameter_df,
        observable_df = observable_df
    )
    if petab.lint_problem(problem):
        raise Exception('PEtab lint failed')

    # Write PEtab problem to disk
    if folder is not None:
        folder = os.path.abspath(folder)
        os.makedirs(folder, exist_ok=True)
        problem.to_files(
            sbml_file=os.path.join(folder, f'{modelname}_model.xml'),
            condition_file=os.path.join(folder, f'{modelname}_conditions.tsv'),
            measurement_file=os.path.join(folder, f'{modelname}_measurements.tsv'),
            parameter_file=os.path.join(folder, f'{modelname}_parameters.tsv'),
            observable_file=os.path.join(folder, f'{modelname}_observables.tsv'),
            yaml_file=os.path.join(folder, f'{modelname}.yaml')
        )
        return os.path.join(folder, f'{modelname}.yaml'), T

    else:
        return problem, T


def simulate_splines(*args, folder=None, keep_temporary=False, **kwargs):
    if folder is not None:
        return _simulate_splines(folder, *args, **kwargs)
    elif keep_temporary:
        folder = tempfile.TemporaryDirectory().name
        return _simulate_splines(folder, *args, **kwargs)
    else:
        with tempfile.TemporaryDirectory() as folder:
            return _simulate_splines(folder, *args, **kwargs)


def _default_initial_values(nsplines: int):
    return np.zeros(nsplines)


def random_suffix(n):
    chars = 'abcdefghijklmnopqrstuvwxyz'
    chars += chars.upper()
    chars += '0123456789'
    chars = np.asarray(list(chars))
    chars = chars[np.random.randint(0, len(chars), n)]
    return ''.join(chars)


def _simulate_splines(folder, splines, params_true, initial_values=None, *, rtol=1e-12, atol=1e-12, simulate_upsample=10, discard_annotations=False, use_adjoint=False, skip_sensitivity=False, **kwargs):
    # Default initial values
    if initial_values is None:
        initial_values = _default_initial_values(len(splines))

    # Create PEtab problem
    path, T = create_petab_problem(
        splines, params_true, initial_values,
        sigma=0.0, folder=folder, **kwargs
    )
    problem = petab.Problem.from_yaml(path)

    # Create and compile AMICI model
    model = import_petab_problem(
        problem,
        discard_annotations=discard_annotations,
        model_output_dir=os.path.join(folder, 'amici_models'),
        model_name='splinetest_' + random_suffix(6)
                                   # to prevent module collisions
    )

    # Set solver options
    solver = model.getSolver()
    solver.setRelativeTolerance(rtol)
    solver.setAbsoluteTolerance(atol)
    if not skip_sensitivity:
        solver.setSensitivityOrder(amici.SensitivityOrder_first)
        if use_adjoint:
            solver.setSensitivityMethod(amici.SensitivityMethod_adjoint)
        else:
            solver.setSensitivityMethod(amici.SensitivityMethod_forward)

    # Compute and set timepoints
    n = max(len(spline.xx) for spline in splines) * simulate_upsample
    tt = np.linspace(0, float(T), n)
    model.setTimepoints(tt)

    # Simulate PEtab problem
    params_str = {p.name : v for p, v in params_true.items()}
    res = simulate_petab(problem, model, solver, params_str)
    llh, sllh, rdatas = res[LLH], res[SLLH], res[RDATAS]
    assert len(rdatas) == 1
    rdata = rdatas[0]

    # Return state/parameter ordering
    state_ids = model.getStateIds()
    param_ids = model.getParameterIds()

    return llh, sllh, rdata, state_ids, param_ids


def check_splines(splines,
                  params_true,
                  initial_values=None,
                  *,
                  discard_annotations=False,
                  use_adjoint=False,
                  skip_sensitivity=False,
                  debug=False,
                  llh_rtol=1e-8,
                  sllh_atol=1e-8,
                  x_rtol=1e-11,
                  x_atol=1e-11,
                  w_rtol=1e-11,
                  w_atol=1e-11,
                  sx_rtol=1e-10,
                  sx_atol=1e-10,
                  **kwargs
    ):

    if isinstance(splines, AbstractSpline):
        splines = [splines]

    # Default initial values
    if initial_values is None:
        initial_values = _default_initial_values(len(splines))

    # Simulate PEtab problem
    llh, sllh, rdata, state_ids, param_ids = simulate_splines(
        splines, params_true, initial_values,
        discard_annotations=discard_annotations,
        skip_sensitivity=skip_sensitivity,
        use_adjoint=use_adjoint,
        **kwargs
    )

    tt = rdata['ts']

    # Sort splines/ics/parameters as in the AMICI model
    splines = [splines[species_to_index(name)] for name in state_ids]
    initial_values = [initial_values[species_to_index(name)] for name in state_ids]
    def param_by_name(id):
        for p in params_true.keys():
            if p.name == id:
                return p
        assert False
    params_sorted = [param_by_name(id) for id in param_ids]

    # Check states
    x_true_sym = sp.Matrix([
        integrate_spline(spline, None, tt, iv)
        for (spline, iv) in zip(splines, initial_values)
    ]).transpose()
    x_true = np.asarray(x_true_sym.subs(params_true), dtype=float)
    if not debug:
        check_results(rdata, 'x', x_true, assert_fun, x_atol, x_rtol)
    elif debug == 'print':
        x_err_abs = abs(rdata['x'] - x_true)
        x_err_rel = np.where(
            x_err_abs == 0,
            0,
            x_err_abs / abs(x_true)
        )
        print("x_err_abs:")
        print(np.squeeze(x_err_abs))
        print("x_err_rel:")
        print(np.squeeze(x_err_rel))

    # Check spline evaluations
    # TODO can we know how the splines are ordered inside w?
    if False and discard_annotations and len(splines) == 1:
        assert rdata['w'].shape[1] == 1
        w_true = np.column_stack([
            evaluate_spline(spline, params_true, tt, dtype=float)
            for spline in splines
        ])
        if not debug:
            check_results(rdata, 'w', w_true, assert_fun, w_atol, w_rtol)
        elif debug == 'print':
            w_err_abs = abs(rdata['w'] - w_true)
            w_err_rel = np.where(
                w_err_abs == 0,
                0,
                w_err_abs / abs(w_true)
            )
            print("w_err_abs:")
            print(np.squeeze(w_err_abs))
            print("w_err_rel:")
            print(np.squeeze(w_err_rel))
    else:
        w_true = None

    # Check sensitivities
    sx_by_state = [
        x_true_sym[:, i].jacobian(params_sorted).subs(params_true)
        for i in range(x_true_sym.shape[1])
    ]
    sx_by_state = [np.asarray(sx, dtype=float) for sx in sx_by_state]
    sx_true = np.concatenate([
        sx[:, :, np.newaxis] for sx in sx_by_state
    ], axis=2)
    if skip_sensitivity or use_adjoint:
        pass
    elif not debug:
        check_results(rdata, 'sx', sx_true, assert_fun, sx_atol, sx_rtol)
    elif debug == 'print':
        sx_err_abs = abs(rdata['sx'] - sx_true)
        sx_err_rel = np.where(
            sx_err_abs == 0,
            0,
            sx_err_abs / abs(sx_true)
        )
        print("sx_err_abs:")
        print(np.squeeze(sx_err_abs))
        print("sx_err_rel:")
        print(np.squeeze(sx_err_rel))

    # Check log-likelihood
    llh_true = - 0.5 * rdata['y'].size * np.log(2*np.pi)
    llh_error_rel = abs(llh - llh_true) / abs(llh_true)
    if not debug:
        assert llh_error_rel <= llh_rtol
    elif debug == 'print':
        print(f'llh_error_rel = {llh_error_rel}')

    # Check log-likelihood sensitivities
    # (should be all zero, since we simulated with the true parameters)
    if not skip_sensitivity:
        if sllh_atol is None:
            sllh_atol = np.finfo(float).eps
        sllh = np.asarray([s for s in sllh.values()])
        sllh_err_abs = abs(sllh).max()
        if not debug:
            assert sllh_err_abs <= sllh_atol
        elif debug == 'print':
            print(f'sllh_err_abs = {sllh_err_abs}')

    if debug:
        return dict(
            rdata=rdata,
            splines=splines,
            initial_values=initial_values,
            params_true=params_true,
            params_sorted=params_sorted,
            x_true=x_true,
            w_true=w_true,
            sx_true=sx_true,
            llh_true=llh_true
        )


def check_splines_full(*args, **kwargs):
    check_splines(*args, **kwargs, discard_annotations=True, use_adjoint=False)
    check_splines(*args, **kwargs, discard_annotations=False, use_adjoint=False)
    check_splines(*args, **kwargs, discard_annotations=False, use_adjoint=True)


def example_spline_1(idx: int = 0, offset: float = 0, scale: float = 1, num_nodes: int = 9, fixed_values=None):
    "A simple spline with no extrapolation."

    yy_true = np.asarray([0.0, 2.0, 5.0, 6.0, 5.0, 4.0, 2.0, 3.0, 4.0, 6.0, 7.0, 7.5, 6.5, 4.0])
    if num_nodes is not None:
        assert 1 < num_nodes <= len(yy_true)
        yy_true = yy_true[:num_nodes]
    yy_true = scale * yy_true + offset

    xx = UniformGrid(0, 25, length=len(yy_true))

    yy = list(sp.symbols(f'y{idx}_0:{len(yy_true)}'))

    if fixed_values is None:
        params = dict(zip(yy, yy_true))
    else:
        params = {}
        for i in range(len(yy_true)):
            if i in fixed_values:
                yy[i] = yy_true[i]
            else:
                params[yy[i]] = yy_true[i]

    spline = CubicHermiteSpline(
        f'y{idx}', amici_time_symbol,
        xx, yy,
        bc=None, extrapolate=None
    )

    tols = dict(llh_rtol=1e-15)

    return spline, params, tols


def example_spline_2(idx: int = 0):
    "A simple spline with periodic bc but no extrapolation."
    yy_true = [0.0, 2.0, 3.0, 4.0, 1.0, -0.5, -1, -1.5, 0.5, 0.0]
    xx = UniformGrid(0, 25, length=len(yy_true))
    yy = list(sp.symbols(f'y{idx}_0:{len(yy_true) - 1}'))
    yy.append(yy[0])
    params = dict(zip(yy, yy_true))
    spline = CubicHermiteSpline(
        f'y{idx}', amici_time_symbol,
        xx, yy,
        bc='periodic', extrapolate=None
    )
    tols = {}
    return spline, params, tols


def example_spline_3(idx: int = 0):
    "A simple spline with extrapolation on the right side."
    yy_true = [0.0, 2.0, 5.0, 6.0, 5.0, 4.0, 2.0, 3.0, 4.0, 6.0]
    xx = UniformGrid(0, 25, length=len(yy_true))
    yy = list(sp.symbols(f'y{idx}_0:{len(yy_true)}'))
    params = dict(zip(yy, yy_true))
    spline = CubicHermiteSpline(
        f'y{idx}', amici_time_symbol,
        xx, yy,
        bc=None, extrapolate=(None, 'constant')
    )
    tols = dict(llh_rtol=1e-15)
    return spline, params, tols


def example_splines_1():
    spline0, params0, tols0 = example_spline_1(0, num_nodes=9, fixed_values=[0, 2])
    spline1, params1, tols1 = example_spline_1(1, num_nodes=14, scale=1.5, offset=5)
    spline2, params2, tols2 = example_spline_1(2, num_nodes=5,  scale=0.5, offset=-5)

    splines = [spline0, spline1, spline2]

    params = dict(params0)
    params.update(params1)
    params.update(params2)

    keys = set().union(tols0.keys(), tols1.keys(), tols2.keys())
    tols = {}
    for key in keys:
        tols[key] = max(
            tols0[key] if key in tols0.keys() else 0.0,
            tols1[key] if key in tols1.keys() else 0.0,
            tols2[key] if key in tols2.keys() else 0.0,
        )

    return splines, params, tols


def test_CubicHermiteSpline():
    spline, params, tols = example_spline_1()
    check_splines_full(spline, params, **tols)

    # same as above, but with some fixed values
    spline, params, tols = example_spline_1(fixed_values=[0, 2])
    check_splines_full(spline, params, **tols)

    spline, params, tols = example_spline_2()
    check_splines_full(spline, params, **tols)


def test_multiple_splines():
    splines, params, tols = example_splines_1()
    check_splines_full(splines, params, **tols)


def test_splines_evaluated_at_formula():
    pass
