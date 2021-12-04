import math
import os
import tempfile
import uuid
from typing import List, Optional, Union

import numpy as np
import pandas as pd
import petab
import sympy as sp

import amici
from amici.gradient_check import check_results
from amici.petab_import import import_petab_problem
from amici.petab_objective import (LLH, RDATAS, SLLH, simulate_petab)
from amici.sbml_utils import (addCompartment, addInflow, addParameter,
                              addRateRule, addSpecies, amici_time_symbol,
                              createSbmlModel)
from amici.splines import AbstractSpline, CubicHermiteSpline, UniformGrid


# callback for amici.gradient_check.check_*
def assert_fun(x):
    assert x


def evaluate_spline(spline: AbstractSpline, params: dict, tt, **kwargs):
    """
    Evaluate the `AbstractSpline` `spline` at timepoints `tt`
    for the parameters given in the dictionary `params`.
    """
    return np.asarray([spline.evaluate(t).subs(params) for t in tt], **kwargs)


def integrate_spline(
        spline: AbstractSpline, params: dict, tt, initial_value=0, **kwargs
        ):
    """
    Integrate the `AbstractSpline` `spline` at timepoints `tt`
    for the parameters given in the dictionary `params`.
    """
    ispline = [initial_value + spline.integrate(0, t) for t in tt]
    if params is not None:
        ispline = [x.subs(params) for x in ispline]
    return np.asarray(ispline, **kwargs)


def create_condition_table():
    """Create a PEtab condition table."""
    condition_df = pd.DataFrame({'conditionId': ['condition1']})
    condition_df.set_index(['conditionId'], inplace=True)
    return condition_df


def create_parameter_table(**columns):
    """Create a PEtab parameter table."""
    if isinstance(columns['parameterId'], str):
        columns['parameterId'] = [columns['parameterId']]
    columns.setdefault('parameterScale', 'lin')
    columns.setdefault('estimate', 1)
    parameter_df = pd.DataFrame(columns)
    parameter_df.set_index(['parameterId'], inplace=True)
    return parameter_df


def create_observable_table(**columns):
    """Create a PEtab observable table."""
    if isinstance(columns['observableId'], str):
        columns['observableId'] = [columns['observableId']]
    columns.setdefault('observableTransformation', 'lin')
    columns.setdefault('noiseDistribution', 'normal')
    observable_df = pd.DataFrame(columns)
    observable_df.set_index(['observableId'], inplace=True)
    return observable_df


def create_measurement_table(**columns):
    """Create a PEtab measurement table."""
    if isinstance(columns['observableId'], str):
        columns['observableId'] = [columns['observableId']]
    columns.setdefault('simulationConditionId', 'condition1')
    return pd.DataFrame(columns)


def species(i):
    return f'z{i}'


def observable(i):
    return f'{species(i)}_obs'


def species_to_index(name):
    assert name[0] == 'z'
    return int(name[1:])


def create_petab_problem(
        splines: List[AbstractSpline],
        params_true,
        initial_values,
        use_reactions=False,
        measure_upsample=6,
        sigma=1.0,
        Textrapolate=0.25,
        folder=None,
        model_name='test_splines',
):
    """
    Given a list of `AbstractSplines`, create a PEtab problem for the system of
    ODEs given by `z_i(t)' = spline[i](t)`.

    :param params_true:
        parameter values used to compute the analytical solution of the ODE
        system in order to fill the PEtab measurement table

    :param initial_values:
        initial values of the state variables

    :param use_reactions:
        whether the ODEs are encoded in the SBML model as reactions (inflows)
        or rate rules

    :param measure_upsample:
        controls the number of time points at which synthetic measurements are
        taken. The interval between subsequent time points is equal to the
        smallest interval between subsequent spline nodes divided by
        `measure_upsample`

    :param sigma:
        standard deviation for additive Normal noise used to corrupt synthetic
        measurements

    :param Textrapolate:
        factor controlling how long after the final spline node the simulation
        should continue in order to test extrapolation methods.

    :param folder:
        if not `None`, save the PEtab problem to this folder

    :param model_name:
        name of the SBML model to be created
    """

    for spline in splines:
        if spline.x != amici_time_symbol:
            raise Exception(
                'the given splines must be evaluated at the simulation time'
            )

    # Create SBML document
    doc, model = createSbmlModel(model_name)
    addCompartment(model, 'compartment')
    for (i, spline) in enumerate(splines):
        spline.addToSbmlModel(model)
        addSpecies(model, species(i), initial_amount=initial_values[i])
        if use_reactions:
            addInflow(model, species(i), splines[i].sbml_id)
        else:
            addRateRule(model, species(i), splines[i].sbml_id)
    for (parId, value) in params_true.items():
        addParameter(model, parId, value=value, constant=True)
    for spline in splines:
        addParameter(model, spline.sbml_id, constant=False)

    # Compute simulation time
    # Must cover all the intervals of definition for the splines,
    # plus something extra for extrapolated or periodic splines
    T = 0
    for spline in splines:
        if spline.extrapolate[0] is None and spline.xx[0] > 0:
            raise Exception(
                'if no left-extrapolation is defined for a spline, '
                'its interval of definition should contain zero'
            )
        if spline.extrapolate[1] is not None:
            f = Textrapolate if spline.extrapolate[
                                    1] != 'periodic' else 1 + Textrapolate
            DT = f * (spline.xx[-1] - spline.xx[0])
        else:
            DT = 0
        T = max(T, spline.xx[-1] + DT)

    # Compute synthetic measurements
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
    # ensure that same parameter order is used for all columns
    _params = list(params_true.items())
    parameter_df = create_parameter_table(
        parameterId=[p.name for (p, v) in _params],
        lowerBound=min(v for (p, v) in _params),
        upperBound=max(v for (p, v) in _params),
        nominalValue=[v for (p, v) in _params],
        estimate=1,
    )
    observable_df = create_observable_table(
        observableId=[observable(i) for i in range(len(splines))],
        observableFormula=[species(i) for i in range(len(splines))],
        noiseFormula=sigma if sigma > 0 else 1.0,
    )
    measurement_df = create_measurement_table(
        observableId=np.concatenate(
            [len(tt_obs) * [observable(i)] for i in range(len(splines))]),
        time=len(splines) * list(tt_obs),
        measurement=np.concatenate(zz_obs),
    )

    # Create and validate PEtab problem
    problem = petab.Problem(
        sbml_document=doc,
        sbml_model=model,
        condition_df=condition_df,
        measurement_df=measurement_df,
        parameter_df=parameter_df,
        observable_df=observable_df,
    )
    if petab.lint_problem(problem):
        raise Exception('PEtab lint failed')

    # Write PEtab problem to disk
    if folder is None:
        return problem, T
    folder = os.path.abspath(folder)
    os.makedirs(folder, exist_ok=True)
    problem.to_files(
        sbml_file=os.path.join(folder, f'{model_name}_model.xml'),
        condition_file=os.path.join(folder, f'{model_name}_conditions.tsv'),
        measurement_file=os.path.join(folder,
                                      f'{model_name}_measurements.tsv'),
        parameter_file=os.path.join(folder, f'{model_name}_parameters.tsv'),
        observable_file=os.path.join(folder,
                                     f'{model_name}_observables.tsv'),
        yaml_file=os.path.join(folder, f'{model_name}.yaml'),
    )
    return os.path.join(folder, f'{model_name}.yaml'), T


def simulate_splines(
        splines,
        params_true,
        initial_values=None,
        *,
        folder: Optional[str] = None,
        keep_temporary: bool = False,
        benchmark: Union[bool, int] = False,
        rtol: float = 1e-12,
        atol: float = 1e-12,
        maxsteps: int = 500_000,
        discard_annotations: bool = False,
        use_adjoint: bool = False,
        skip_sensitivity: bool = False,
        **kwargs
):
    """
    Create a PEtab problem using `create_petab_problem` and simulate it with
    AMICI.

    :param splines:
        passed to `create_petab_problem`

    :param params_true:
        passed to `create_petab_problem`

    :param initial_values:
        passed to `create_petab_problem`

    :param folder:
        working directory (a temporary one is used if not specified)

    :param keep_temporary:
        whether to keep or delete temporary working directories on exit

    :param benchmark:
        instead of returning the simulation data, run the simulation
        `benchmark` times (defaults to `50` if `benchmark` is `True`)
        and return execution times

    :param rtol:
        relative tolerance for AMICI solver

    :param atol:
        absolute tolerance for AMICI solver

    :param maxsteps:
        maximum number of steps for AMICI solver

    :param discard_annotations:
        whether to discard spline annotations,
        forcing AMICI to read the spline as a piecewise assignment rule

    :param use_adjoint:
        whether to use adjoint sensitivity computation

    :param skip_sensitivity:
        whether to skip sensitivity computation

    :param kwargs:
        passed to `create_petab_problem`
    """

    # If no working directory is given, create a temporary one
    if folder is None:
        if keep_temporary:
            folder = tempfile.TemporaryDirectory().name
            print(f"temporary folder is {folder}")
        else:
            with tempfile.TemporaryDirectory() as folder:
                return simulate_splines(
                    splines, params_true, initial_values, folder=folder,
                    benchmark=benchmark, rtol=rtol, atol=atol,
                    maxsteps=maxsteps, discard_annotations=discard_annotations,
                    use_adjoint=use_adjoint, skip_sensitivity=skip_sensitivity,
                    **kwargs
                )

    # Default initial values
    if initial_values is None:
        initial_values = np.zeros(len(splines))

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
        model_name='splinetest_' + uuid.uuid1().hex
        # to prevent module collisions
    )

    # Set solver options
    solver = model.getSolver()
    solver.setRelativeTolerance(rtol)
    solver.setAbsoluteTolerance(atol)
    solver.setMaxSteps(maxsteps)
    if not skip_sensitivity:
        solver.setSensitivityOrder(amici.SensitivityOrder_first)
        if use_adjoint:
            solver.setSensitivityMethod(amici.SensitivityMethod_adjoint)
        else:
            solver.setSensitivityMethod(amici.SensitivityMethod_forward)

    # Compute and set timepoints
    # NB not working, will always be equal to the observation times
    # n = max(len(spline.xx) for spline in splines) * simulate_upsample
    # tt = np.linspace(0, float(T), n)
    # model.setTimepoints(tt)

    # Create dictionary for parameter values
    params_str = {p.name: v for (p, v) in params_true.items()}

    if benchmark is False:
        # Simulate PEtab problem
        res = simulate_petab(problem, model, solver, params_str)
        assert SLLH not in res.keys()
        llh, rdatas = res[LLH], res[RDATAS]
        assert len(rdatas) == 1
        llh = float(llh)
        rdata = rdatas[0]
        assert SLLH in rdata.keys()
        sllh = rdata[SLLH]

        # Return state/parameter ordering
        state_ids = model.getStateIds()
        param_ids = model.getParameterIds()

        return initial_values, llh, sllh, rdata, state_ids, param_ids

    else:
        if benchmark is True:
            benchmark = 50
        import time
        runtimes = []
        for _ in range(int(benchmark)):
            t0 = time.perf_counter()
            simulate_petab(problem, model, solver, params_str)
            t_elapsed = time.perf_counter() - t0
            runtimes.append(t_elapsed)
        return dict(
            runtimes=runtimes,
            mean=np.mean(runtimes),
            median=np.median(runtimes),
            min=min(runtimes),
            max=max(runtimes),
        )


def check_splines(
        splines,
        params_true,
        initial_values=None,
        *,
        assert_fun,
        discard_annotations: bool = False,
        use_adjoint: bool = False,
        skip_sensitivity: bool = False,
        debug: bool = False,
        llh_rtol: float = 1e-8,
        sllh_atol: float = 1e-8,
        x_rtol: float = 1e-11,
        x_atol: float = 1e-11,
        w_rtol: float = 1e-11,
        w_atol: float = 1e-11,
        sx_rtol: float = 1e-10,
        sx_atol: float = 1e-10,
        **kwargs
):
    """
    Create a PEtab problem using `create_petab_problem`,
    simulate it with `simulate_splines`
    and check it against the analytical solution.

    :param splines:
        passed to `simulate_splines`

    :param params_true:
        passed to `simulate_splines`

    :param initial_values:
        passed to `simulate_splines`

    :param assert_fun:
        function to use for checks

    :param discard_annotations:
        whether to discard spline annotations,
        forcing AMICI to read the spline as a piecewise assignment rule

    :param use_adjoint:
        whether to use adjoint sensitivity computation

    :param skip_sensitivity:
        whether to skip sensitivity computation

    :param debug:
        if not `False`, do not check and return results and ground truth
        instead.
        If equal to `'print'`, in addition to the above print error values.

    :param kwargs:
        passed to `simulate_splines`
    """

    if isinstance(splines, AbstractSpline):
        splines = [splines]

    # Simulate PEtab problem
    initial_values, llh, sllh, rdata, state_ids, param_ids = simulate_splines(
        splines, params_true, initial_values,
        discard_annotations=discard_annotations,
        skip_sensitivity=skip_sensitivity,
        use_adjoint=use_adjoint,
        **kwargs
    )

    tt = rdata['ts']

    # Sort splines/ics/parameters as in the AMICI model
    splines = [splines[species_to_index(name)] for name in state_ids]
    initial_values = [initial_values[species_to_index(name)] for name in
                      state_ids]

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
    llh_true = - 0.5 * rdata['y'].size * np.log(2 * np.pi)
    llh_error_rel = abs(llh - llh_true) / abs(llh_true)
    if not debug:
        assert_fun(llh_error_rel <= llh_rtol)
    elif debug == 'print':
        print(f'llh_error_rel = {llh_error_rel}')

    # Check log-likelihood sensitivities
    # (should be all zero, since we simulated with the true parameters)
    if not skip_sensitivity:
        if sllh_atol is None:
            sllh_atol = np.finfo(float).eps
        sllh_err_abs = abs(sllh).max()
        if not debug:
            assert_fun(sllh_err_abs <= sllh_atol)
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


def check_splines_full(
        splines, params, tols, *args, check_piecewise=True, **kwargs
        ):
    """
    Check example PEtab problem with `check_splines`
    both using adjoint and forward sensitivities
    and also in the case in which the splines are read as piecewise functions.
    """
    if isinstance(tols, dict):
        tols1 = tols2 = tols3 = tols
    else:
        tols1, tols2, tols3 = tols

    if isinstance(splines, AbstractSpline):
        contains_periodic = (splines.extrapolate == ('periodic', 'periodic'))
    elif any(spline.extrapolate == ('periodic', 'periodic') for spline in
             splines):
        contains_periodic = True
    else:
        contains_periodic = False

    if check_piecewise and not contains_periodic:
        check_splines(splines, params, *args, **kwargs, **tols1,
                      discard_annotations=True, use_adjoint=False)
    check_splines(splines, params, *args, **kwargs, **tols2,
                  discard_annotations=False, use_adjoint=False)
    check_splines(splines, params, *args, **kwargs, **tols3,
                  discard_annotations=False, use_adjoint=True)


def example_spline_1(
        idx: int = 0,
        offset: float = 0,
        scale: float = 1,
        num_nodes: int = 9,
        fixed_values=None,  # a list of indices or 'all'
):
    """A simple spline with no extrapolation."""

    yy_true = np.asarray(
        [0.0, 2.0, 5.0, 6.0, 5.0, 4.0, 2.0, 3.0, 4.0, 6.0, 7.0, 7.5, 6.5, 4.0])
    if num_nodes is not None:
        assert 1 < num_nodes <= len(yy_true)
        yy_true = yy_true[:num_nodes]
    yy_true = scale * yy_true + offset

    xx = UniformGrid(0, 25, length=len(yy_true))

    yy = list(sp.symbols(f'y{idx}_0:{len(yy_true)}'))

    if fixed_values is None:
        params = dict(zip(yy, yy_true))
    elif fixed_values == 'all':
        params = {}
        for i in range(len(yy_true)):
            yy[i] = yy_true[i]
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
    """A simple spline with periodic bc but no extrapolation."""
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
    tols = (
        dict(llh_rtol=1e-15),
        dict(llh_rtol=1e-15),
        dict(llh_rtol=1e-15, sllh_atol=5e-8, x_rtol=1e-10, x_atol=5e-10)
    )
    return spline, params, tols


def example_spline_3(idx: int = 0):
    """A simple spline with extrapolation on the right side."""
    yy_true = [0.0, 2.0, 5.0, 6.0, 5.0, 4.0, 2.0, 3.0, 4.0, 6.0]
    xx = UniformGrid(0, 25, length=len(yy_true))
    yy = list(sp.symbols(f'y{idx}_0:{len(yy_true)}'))
    params = dict(zip(yy, yy_true))
    spline = CubicHermiteSpline(
        f'y{idx}', amici_time_symbol,
        xx, yy,
        bc=(None, 'zeroderivative'), extrapolate=(None, 'constant')
    )
    tols = {}
    return spline, params, tols


def example_spline_4(idx: int = 0):
    """A simple spline with periodic extrapolation."""
    yy_true = [0.0, 2.0, 3.0, 4.0, 1.0, -0.5, -1, -1.5, 0.5, 0.0]
    xx = UniformGrid(0, 25, length=len(yy_true))
    yy = list(sp.symbols(f'y{idx}_0:{len(yy_true) - 1}'))
    yy.append(yy[0])
    params = dict(zip(yy, yy_true))
    spline = CubicHermiteSpline(
        f'y{idx}', amici_time_symbol,
        xx, yy,
        bc='periodic', extrapolate='periodic'
    )
    tols = dict(llh_rtol=1e-15, sllh_atol=5e-7, x_rtol=5e-10, x_atol=7.5e-10,
                sx_rtol=5e-10, sx_atol=5e-10)
    return spline, params, tols


def example_splines_1():
    spline0, params0, tols0 = example_spline_1(0, num_nodes=9,
                                               fixed_values=[0, 2])
    spline1, params1, tols1 = example_spline_1(1, num_nodes=14, scale=1.5,
                                               offset=5)
    spline2, params2, tols2 = example_spline_1(2, num_nodes=5, scale=0.5,
                                               offset=-5)

    splines = [spline0, spline1, spline2]

    params = dict(params0)
    params.update(params1)
    params.update(params2)

    keys = set().union(tols0.keys(), tols1.keys(), tols2.keys())
    tols = {key: max(
        tols0[key] if key in tols0.keys() else 0.0,
        tols1[key] if key in tols1.keys() else 0.0,
        tols2[key] if key in tols2.keys() else 0.0,
    ) for key in keys}
    tols['llh_rtol'] = 1e-14
    tols['sllh_atol'] = 5e-8
    tols['sx_rtol'] = 1e-9
    tols['sx_atol'] = 1e-10
    tols['x_rtol'] = 1e-10
    tols['x_atol'] = 1e-9

    return splines, params, tols


def test_CubicHermiteSpline(**kwargs):
    spline, params, tols = example_spline_1()
    check_splines_full(spline, params, tols, assert_fun=assert_fun, **kwargs)

    # same as above, but with some fixed values
    spline, params, tols = example_spline_1(fixed_values=[0, 2])
    check_splines_full(spline, params, tols, assert_fun=assert_fun, **kwargs)

    # same as above, but with all values fixed
    spline, params, tols = example_spline_1(fixed_values='all')
    check_splines_full(spline, params, tols, assert_fun=assert_fun, **kwargs)

    spline, params, tols = example_spline_2()
    check_splines_full(spline, params, tols, assert_fun=assert_fun, **kwargs)


def test_multiple_splines(**kwargs):
    splines, params, tols = example_splines_1()
    check_splines_full(splines, params, tols, assert_fun=assert_fun, **kwargs)


def test_splines_evaluated_at_formula():
    raise NotImplementedError("Implement me!")
