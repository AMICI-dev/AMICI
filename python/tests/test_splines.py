import math
import os
import uuid
from typing import List, Optional, Union, Dict, Sequence

import amici
import petab
import pytest
import sympy as sp
from amici.gradient_check import _check_results
from amici.petab_import import import_petab_problem
from amici.petab_objective import (LLH, RDATAS, SLLH, EDATAS, simulate_petab)
from amici.sbml_utils import (add_compartment, add_inflow, add_parameter,
                              add_rate_rule, add_species, amici_time_symbol,
                              create_sbml_model)
from amici.splines import AbstractSpline, CubicHermiteSpline, UniformGrid
from amici.testing import TemporaryDirectoryWinSafe as TemporaryDirectory

import numpy as np
import pandas as pd


def evaluate_spline(spline: AbstractSpline, params: dict, tt, **kwargs):
    """
    Evaluate the `AbstractSpline` `spline` at timepoints `tt`
    for the parameters given in the dictionary `params`.
    """
    return np.asarray([spline.evaluate(t).subs(params) for t in tt], **kwargs)


def integrate_spline(
        spline: AbstractSpline, params: Union[Dict, None], tt, initial_value=0,
        **kwargs
):
    """
    Integrate the `AbstractSpline` `spline` at timepoints `tt`
    for the parameters given in the dictionary `params`.
    """
    ispline = [initial_value + spline.integrate(0, t) for t in tt]
    if params is not None:
        ispline = [x.subs(params) for x in ispline]
    return np.asarray(ispline, **kwargs)


def create_condition_table() -> pd.DataFrame:
    """Create a PEtab condition table."""
    condition_df = pd.DataFrame({'conditionId': ['condition1']})
    condition_df.set_index(['conditionId'], inplace=True)
    return condition_df


def create_parameter_table(**columns) -> pd.DataFrame:
    """Create a PEtab parameter table."""
    if isinstance(columns['parameterId'], str):
        columns['parameterId'] = [columns['parameterId']]
    columns.setdefault('parameterScale', 'lin')
    columns.setdefault('estimate', 1)
    parameter_df = pd.DataFrame(columns)
    parameter_df.set_index(['parameterId'], inplace=True)
    return parameter_df


def create_observable_table(**columns) -> pd.DataFrame:
    """Create a PEtab observable table."""
    if isinstance(columns['observableId'], str):
        columns['observableId'] = [columns['observableId']]
    columns.setdefault('observableTransformation', 'lin')
    columns.setdefault('noiseDistribution', 'normal')
    observable_df = pd.DataFrame(columns)
    observable_df.set_index(['observableId'], inplace=True)
    return observable_df


def create_measurement_table(**columns) -> pd.DataFrame:
    """Create a PEtab measurement table."""
    if isinstance(columns['observableId'], str):
        columns['observableId'] = [columns['observableId']]
    columns.setdefault('simulationConditionId', 'condition1')
    return pd.DataFrame(columns)


def species(i) -> str:
    return f'z{i}'


def observable(i) -> str:
    return f'{species(i)}_obs'


def species_to_index(name) -> int:
    assert name[0] == 'z'
    return int(name[1:])


def create_petab_problem(
        splines: List[AbstractSpline],
        params_true,
        initial_values=None,
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
            raise ValueError(
                'the given splines must be evaluated at the simulation time'
            )

    if initial_values is None:
        initial_values = np.zeros(len(splines))

    # Create SBML document
    doc, model = create_sbml_model(model_name)
    add_compartment(model, 'compartment')
    for (i, spline) in enumerate(splines):
        spline.add_to_sbml_model(model)
        add_species(model, species(i), initial_amount=initial_values[i])
        if use_reactions:
            add_inflow(model, species(i), splines[i].sbml_id)
        else:
            add_rate_rule(model, species(i), splines[i].sbml_id)
    for (parId, value) in params_true.items():
        add_parameter(model, parId, value=value, constant=True)
    for spline in splines:
        add_parameter(model, spline.sbml_id, constant=False)

    # Compute simulation time
    # Must cover all the intervals of definition for the splines,
    # plus something extra for extrapolated or periodic splines
    T = 0
    for spline in splines:
        if spline.extrapolate[0] is None and spline.xx[0] > 0:
            raise ValueError(
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
        lowerBound=min(v for (p, v) in _params) if _params else [],
        upperBound=max(v for (p, v) in _params) if _params else [],
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
        raise RuntimeError('PEtab lint failed')

    # Write PEtab problem to disk
    if folder is None:
        return problem, initial_values, T
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
    return os.path.join(folder, f'{model_name}.yaml'), initial_values, T


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
        petab_problem=None,
        amici_model=None,
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

    :param petab_problem:
        PEtab problem (if already created)

    :param amici_model:
        AMICI model (if already created)

    :param kwargs:
        passed to `create_petab_problem`
    """
    # If no working directory is given, create a temporary one
    if folder is None:
        if keep_temporary:
            folder = TemporaryDirectory().name
            print(f"temporary folder is {folder}")
        else:
            with TemporaryDirectory() as folder:
                return simulate_splines(
                    splines, params_true, initial_values, folder=folder,
                    benchmark=benchmark, rtol=rtol, atol=atol,
                    maxsteps=maxsteps, discard_annotations=discard_annotations,
                    use_adjoint=use_adjoint, skip_sensitivity=skip_sensitivity,
                    petab_problem=petab_problem, amici_model=amici_model,
                    **kwargs
                )

    if petab_problem is None and amici_model is not None:
        raise ValueError(
            "if amici_model is given, petab_problem must be given too"
        )

    if petab_problem is not None and initial_values is None:
        raise ValueError(
            "if petab_problem is given, initial_values must be given too"
        )

    if petab_problem is None:
        # Create PEtab problem
        path, initial_values, T = create_petab_problem(
            splines, params_true, initial_values,
            sigma=0.0, folder=folder, **kwargs
        )
        petab_problem = petab.Problem.from_yaml(path)

    if amici_model is None:
        # Create and compile AMICI model
        amici_model = import_petab_problem(
            petab_problem,
            discard_sbml_annotations=discard_annotations,
            model_output_dir=os.path.join(folder, 'amici_models'),
            # to prevent module collisions
            model_name='splinetest_' + uuid.uuid1().hex
        )

    # Set solver options
    solver = amici_model.getSolver()
    solver.setRelativeTolerance(rtol)
    solver.setAbsoluteTolerance(atol)
    solver.setMaxSteps(maxsteps)
    if not skip_sensitivity:
        solver.setSensitivityOrder(amici.SensitivityOrder.first)
        if use_adjoint:
            solver.setSensitivityMethod(amici.SensitivityMethod.adjoint)
        else:
            solver.setSensitivityMethod(amici.SensitivityMethod.forward)

    # Compute and set timepoints
    # NB not working, will always be equal to the observation times
    # n = max(len(spline.xx) for spline in splines) * simulate_upsample
    # tt = np.linspace(0, float(T), n)
    # model.setTimepoints(tt)

    # Create dictionary for parameter values
    params_str = {p.name: v for (p, v) in params_true.items()}

    if benchmark is False:
        # Simulate PEtab problem
        res = simulate_petab(petab_problem, amici_model, solver, params_str)
        assert SLLH not in res.keys()
        llh, rdatas, edatas = res[LLH], res[RDATAS], res[EDATAS]
        assert len(rdatas) == 1
        llh = float(llh)
        rdata = rdatas[0]
        assert SLLH in rdata.keys()
        sllh = rdata[SLLH]
        assert len(edatas) == 1
        edata = edatas[0]

        # Return state/parameter ordering
        state_ids = amici_model.getStateIds()
        param_ids = amici_model.getParameterIds()

        return initial_values, petab_problem, amici_model, solver, llh, sllh, rdata, edata, state_ids, param_ids

    if benchmark is True:
        benchmark = 50
    import time
    runtimes = []
    for _ in range(int(benchmark)):
        t0 = time.perf_counter()
        simulate_petab(petab_problem, amici_model, solver, params_str)
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
        discard_annotations: bool = False,
        use_adjoint: bool = False,
        skip_sensitivity: bool = False,
        debug: Union[bool, str] = False,
        parameter_lists: Optional[Sequence[Sequence[int]]] = None,
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

    :param parameter_lists:
        Set AMICI parameter list to these values,
        in order to check that partial sensitivity computation works.

    :param kwargs:
        passed to `simulate_splines`
    """

    if isinstance(splines, AbstractSpline):
        splines = [splines]

    # Simulate PEtab problem
    initial_values, petab_problem, amici_model, amici_solver, llh, sllh, rdata, edata, state_ids, param_ids = simulate_splines(
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
        _check_results(rdata, 'x', x_true, atol=x_atol, rtol=x_rtol)
    elif debug == 'print':
        x_err_abs = abs(rdata['x'] - x_true)
        x_err_rel = np.where(
            x_err_abs == 0,
            0,
            x_err_abs / abs(x_true)
        )
        print(f'x_atol={x_atol} x_rtol={x_rtol}')
        print("x_err_abs:")
        print(np.squeeze(x_err_abs))
        print("x_err_abs (maximum):")
        print(x_err_abs.max())
        print("x_err_rel:")
        print(np.squeeze(x_err_rel))
        print("x_err_rel (maximum):")
        print(x_err_rel.max())

    # Check spline evaluations
    # TODO can we know how the splines are ordered inside w?
    if False and discard_annotations and len(splines) == 1:
        assert rdata['w'].shape[1] == 1
        w_true = np.column_stack([
            evaluate_spline(spline, params_true, tt, dtype=float)
            for spline in splines
        ])
        if not debug:
            _check_results(
                rdata, 'w', w_true,
                atol=w_atol, rtol=w_rtol,
            )
        elif debug == 'print':
            w_err_abs = abs(rdata['w'] - w_true)
            w_err_rel = np.where(
                w_err_abs == 0,
                0,
                w_err_abs / abs(w_true)
            )
            print(f'w_atol={w_atol} w_rtol={w_rtol}')
            print("w_err_abs:")
            print(np.squeeze(w_err_abs))
            print("w_err_abs (maximum):")
            print(w_err_abs.max())
            print("w_err_rel:")
            print(np.squeeze(w_err_rel))
            print("w_err_rel (maximum):")
            print(w_err_rel.max())
    else:
        w_true = None

    # Check sensitivities
    if params_sorted:
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
            _check_results(
                rdata, 'sx', sx_true,
                atol=sx_atol, rtol=sx_rtol,
            )
        elif debug == 'print':
            sx_err_abs = abs(rdata['sx'] - sx_true)
            sx_err_rel = np.where(
                sx_err_abs == 0,
                0,
                sx_err_abs / abs(sx_true)
            )
            print(f'sx_atol={sx_atol} sx_rtol={sx_rtol}')
            print("sx_err_abs:")
            print(np.squeeze(sx_err_abs))
            print("sx_err_abs (maximum):")
            print(sx_err_abs.max())
            print("sx_err_rel:")
            print(np.squeeze(sx_err_rel))
            print("sx_err_rel (maximum):")
            print(sx_err_rel.max())
    else:
        assert rdata['sx'] is None

    # Check log-likelihood
    llh_true = - 0.5 * rdata['y'].size * np.log(2 * np.pi)
    llh_error_rel = abs(llh - llh_true) / abs(llh_true)
    if (llh_error_rel > llh_rtol and debug is not True) or debug == 'print':
        print(f'llh_rtol={llh_rtol}')
        print(f'llh_error_rel = {llh_error_rel}')
    if not debug:
        assert llh_error_rel <= llh_rtol

    # Check log-likelihood sensitivities
    # (should be all zero, since we simulated with the true parameters)
    if params_sorted:
        if not skip_sensitivity:
            if sllh_atol is None:
                sllh_atol = np.finfo(float).eps
            sllh_err_abs = abs(sllh).max()
            if (sllh_err_abs > sllh_atol and debug is not True) or debug == 'print':
                print(f'sllh_atol={sllh_atol}')
                print(f'sllh_err_abs = {sllh_err_abs}')
            if not debug:
                assert sllh_err_abs <= sllh_atol
    else:
        assert sllh is None

    # Try different parameter lists
    if not skip_sensitivity and (not use_adjoint) and parameter_lists is not None:
        for plist in parameter_lists:
            amici_model.setParameterList(plist)
            amici_model.setTimepoints(rdata.t)
            rdata_partial = amici.runAmiciSimulation(amici_model, amici_solver)
            assert np.isclose(rdata.sx[:, plist, :], rdata_partial.sx).all()

    if debug:
        return dict(
            splines=splines,
            initial_values=initial_values,
            petab_problem=petab_problem,
            amici_model=amici_model,
            rdata=rdata,
            params_true=params_true,
            params_sorted=params_sorted,
            x_true=x_true,
            w_true=w_true,
            sx_true=sx_true,
            llh_true=llh_true,
        )
    else:
        return dict(
            initial_values=initial_values,
            petab_problem=petab_problem,
            amici_model=amici_model, 
        )


def check_splines_full(
    splines, params, tols, *args,
    check_piecewise=True,
    check_forward=True,
    check_adjoint=True,
    folder=None,
    **kwargs
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
        splines = [splines]

    contains_periodic = any(
        spline.extrapolate == ('periodic', 'periodic')
        for spline in splines
    )

    with TemporaryDirectory() as folder:
        # Amortize creation of PEtab and AMICI objects
        initial_values = None
        petab_problem = None
        amici_model = None

        if check_piecewise and not contains_periodic:
            results = check_splines(
                splines, params, *args, **kwargs, **tols1,
                folder=folder,
                discard_annotations=True,
                use_adjoint=False,
            )
            initial_values = results["initial_values"]
            petab_problem = results["petab_problem"]  

        if check_forward:
            results = check_splines(
                splines, params, *args, **kwargs, **tols2,
                initial_values=initial_values,
                petab_problem=petab_problem,
                use_adjoint=False,
            )
            initial_values = results["initial_values"]
            petab_problem = results["petab_problem"]
            amici_model = results["amici_model"]

        if check_adjoint:
            check_splines(
                splines, params, *args, **kwargs, **tols3,
                initial_values=initial_values, 
                petab_problem=petab_problem,
                amici_model=amici_model,
                use_adjoint=True,
            )


def example_spline_1(
        idx: int = 0,
        offset: float = 0,
        scale: float = 1,
        num_nodes: int = 9,
        fixed_values=None,  # a list of indices or 'all'
        extrapolate=None,
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
        bc=None, extrapolate=extrapolate
    )

    tols = (
        dict(llh_rtol=1e-15),
        dict(llh_rtol=1e-15),
        dict(llh_rtol=1e-15, sllh_atol=5e-8),
    )

    return spline, params, tols


def example_spline_2(idx: int = 0):
    """A simple periodic spline."""
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


def test_spline_piecewise(**kwargs):
    spline, params, tols = example_spline_1()
    check_splines_full(spline, params, tols, **kwargs)


def test_splines(**kwargs):
    spline0, params0, tols0 = example_spline_1(
        0, num_nodes=9, fixed_values=[0, 2], extrapolate='linear'
    )
    spline1, params1, tols1 = example_spline_1(
        1, num_nodes=14, scale=1.5, offset=5, extrapolate='linear'
    )
    spline2, params2, tols2 = example_spline_1(
        2, num_nodes=5, scale=0.5, offset=-5, extrapolate='linear'
    )
    spline3, params3, tols3 = example_spline_1(
        3, fixed_values='all', extrapolate='linear'
    )
    spline4, params4, tols4 = example_spline_2(4)
    spline5, params5, tols5 = example_spline_3(5)

    splines = [spline0, spline1, spline2, spline3, spline4, spline5]

    params = dict(params0)
    params.update(params1)
    params.update(params2)
    params.update(params3)
    params.update(params4)
    params.update(params5)

    if isinstance(tols0, dict):
        tols0 = (tols0, tols0, tols0)
    if isinstance(tols1, dict):
        tols1 = (tols1, tols1, tols1)
    if isinstance(tols2, dict):
        tols2 = (tols2, tols2, tols2)
    if isinstance(tols3, dict):
        tols3 = (tols3, tols3, tols3)
    if isinstance(tols4, dict):
        tols4 = (tols4, tols4, tols4)
    if isinstance(tols5, dict):
        tols5 = (tols5, tols5, tols5)

    tols = []
    for (t0, t1, t2, t3, t4, t5) in zip(tols0, tols1, tols2, tols3, tols4, tols5):
        keys = set().union(t0.keys(), t1.keys(), t2.keys(), t3.keys(), t4.keys(), t5.keys())
        t = {
            key: max(
                t0.get(key, 0.0),
                t1.get(key, 0.0),
                t2.get(key, 0.0),
                t3.get(key, 0.0),
                t4.get(key, 0.0),
                t5.get(key, 0.0),
            ) for key in keys
        }
        tols.append(t)

    tols[1]['x_rtol']   = max(1e-9, tols[1].get('x_rtol', -np.inf))
    tols[1]['x_atol']   = max(5e-9,  tols[1].get('x_atol', -np.inf))
    tols[1]['sx_rtol']  = max(1e-5,  tols[1].get('llh_rtol', -np.inf))
    tols[1]['sx_atol']  = max(5e-9, tols[1].get('sx_atol', -np.inf))
    tols[1]['llh_rtol'] = max(5e-14, tols[1].get('llh_rtol', -np.inf))
    tols[1]['sllh_atol'] = max(5e-5,  tols[1].get('sllh_atol', -np.inf))

    tols[2]['x_rtol']    = max(5e-10, tols[2].get('x_rtol', -np.inf))
    tols[2]['x_atol']    = max(1e-8,  tols[2].get('x_atol', -np.inf))
    tols[2]['llh_rtol']  = max(5e-14, tols[2].get('llh_rtol', -np.inf))
    tols[2]['sllh_atol'] = max(5e-5,  tols[2].get('sllh_atol', -np.inf))

    check_splines_full(splines, params, tols, **kwargs)


def test_splines_plist():
    # Dummy spline #1
    xx = UniformGrid(0, 5, length=3)
    yy = np.asarray([0.0, 1.0, 0.5])
    spline1 = CubicHermiteSpline(
        f'y1', amici_time_symbol,
        xx, yy,
        bc='auto', extrapolate=(None, 'constant'),
    )
    # Dummy spline #2
    xx = UniformGrid(0, 5, length=4)
    yy = np.asarray([0.0, 0.5, -0.5, 0.5])
    spline2 = CubicHermiteSpline(
        f'y2', amici_time_symbol,
        xx, yy,
        bc='auto', extrapolate=(None, 'constant'),
    )
    # Real spline #3
    xx = UniformGrid(0, 5, length=6)
    p1, p2, p3, p4, p5 = sp.symbols('p1 p2 p3 p4 p5')
    yy = np.asarray([p1 + p2, p2 * p3, p4, sp.cos(p1 + p3), p4 * sp.log(p1), p3])
    dd = np.asarray([-0.75, -0.875, p5, 0.125, 1.15057181, 0.0])
    params = {
        p1: 1.0, p2: 0.5, p3: 1.5, p4: -0.25, p5: -0.5
    }
    # print([y.subs(params).evalf() for y in yy])
    spline3 = CubicHermiteSpline(
        f'y3', amici_time_symbol,
        xx, yy, dd,
        bc='auto', extrapolate=(None, 'constant'),
    )
    # Dummy spline 4
    xx = UniformGrid(0, 5, length=3)
    yy = np.asarray([0.0, -0.5, 0.5])
    spline4 = CubicHermiteSpline(
        f'y4', amici_time_symbol,
        xx, yy,
        bc='auto', extrapolate=(None, 'constant'),
    )
    tols = dict(
        x_rtol=1e-6,
        x_atol=1e-11,
        sx_rtol=1e-6,
        sx_atol=5e-11,
        llh_rtol=1e-14,
        sllh_atol=5e-9,
    )
    check_splines_full(
        [spline1, spline2, spline3, spline4], params, tols,
        check_piecewise=False,
        check_forward=False,
        check_adjoint=True, # plist cannot be checked, but complex parameter dependence can
        parameter_lists=[[0, 1, 4], [2, 3]],
    )
    # Debug
    # # res = check_splines(
    #     [spline1, spline2, spline3, spline4], params,
    #     use_adjoint=False,
    #     parameter_lists=[[0, 1, 4], [2, 3]],
    #     #folder='debug',
    #     #debug='print',
    # )
