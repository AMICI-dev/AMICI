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
from amici.splines import CubicHermiteSpline, UniformGrid


# @contextmanager
# def cd(path):
#     cwd = os.getcwd()
#     assert os.path.isabs(cwd)
#     os.chdir(path)
#     yield os.path.abspath(path)
#     os.chdir(pwd)


def assert_fun(x):
    assert x


def integrate_spline(spline, params, tt):
    return np.asarray([spline.integrate(0, t).subs(params) for t in tt])


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
    zz_true = [initial_values[i] + integrate_spline(splines[i], params_true, tt_obs) for i in range(len(splines))]
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


def _simulate_splines(folder, splines, params_true, initial_values=None, *, rtol=1e-12, atol=1e-12, simulate_upsample=10, discard_annotations=False, **kwargs):
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
        #discard_annotations=discard_annotations,
        model_output_dir=os.path.join(folder, 'amici_models')
    )

    # Set solver options
    solver = model.getSolver()
    solver.setRelativeTolerance(rtol)
    solver.setAbsoluteTolerance(atol)
    solver.setSensitivityOrder(amici.SensitivityOrder_first)
    solver.setSensitivityMethod(amici.SensitivityMethod_forward)

    # Compute and set timepoints
    n = max(len(spline.xx) for spline in splines) * simulate_upsample
    tt = np.linspace(0, float(T), n)
    model.setTimepoints(tt)

    # Simulate PEtab problem
    params_str = {p.name : v for p, v in params_true.items()}
    res = simulate_petab(problem, model, solver, params_str)
    llh, sslh, rdatas = res[LLH], res[SLLH], res[RDATAS]
    assert len(rdatas) == 1
    rdata = rdatas[0]

    return llh, sslh, rdata


# def check_splines(splines, params_true, initial_values=None, *, discard_annotations=False, sim_atol=1e-12, sim_rtol=1e-12, check_atol=1e-12, check_rtol=1e-12, sensi_atol=1e-12, fd_atol=1e-12, fd_rtol=1e-12, fd_epsilon=1e-4, **kwargs):
#
#
#         # Set time points for detailed grid
#
#
#
#         # Check that simulation result is correct
#         tt = rdata['ts']
#         for i in range(len(splines)):
#             zz = initial_values[i] + integrate_spline(splines[i], params_true, tt)
#             check_results(rdata, species(i), zz, assert_fun, check_atol, check_rtol)
#
#         # Check that likelihood sensitivities are small near the true parameters
#         # NB they are not zero, not even in exact arithmetic, because the data is finite
#         # TODO check exactly if sigma=zero
#         check_close(sslh, np.zeros_like(), assert_fun, sensi_atol, 0.0, 'sensitivities')
#
#         # Check state/observable sensitivities
#
#         # Check derivatives in a point nearby the true values
#         #params_perturbed = ...
#         edatas = create_parameterized_edatas(model, problem, params_perturbed)
#         assert len(edatas) == 1
#         edata = edatas[0]
#         check_derivatives(model, solver, edata, assert_fun, fd_atol, fd_rtol, fd_epsilon)

#     spline_dt = sp.nsimplify(spline_dt)
#     tmax = (len(yy_true) - 1)*spline_dt
#     xx = UniformGrid(0, tmax, spline_dt)
#     yy = list(sp.symbols(f'y0:{len(yy_true)}'))
#
#     spline = CubicHermiteSpline(
#         'y', amici_time_symbol, xx, yy,
#         bc=bc, extrapolate=extrapolate
#     )
#
# if __name__ == "__main__":
#     import sys
#     folder = sys.argv[1] if len(sys.argv) > 1 else '.'
#     yy_true = [0.0, 2.0, 3.0, 4.0, 1.0, -0.5, -1, -1.5, 0.5, 0.0]
#     create_spline_test_petab(yy_true, folder=folder)
