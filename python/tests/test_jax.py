from pathlib import Path

import amici
import pytest

pytest.importorskip("pysb")
pytest.importorskip("jax")
import diffrax
import jax
import jax.numpy as jnp
import jax.random as jr
import numpy as np
import optimistix
import sympy as sp
from amici import MeasurementChannel as MC
from amici import import_model_module
from amici.exporters.jax.jaxcodeprinter import (
    AmiciJaxCodePrinter,
    _jnp_array_str,
)
from amici.importers.petab.v1 import import_petab_problem
from amici.importers.pysb import pysb2amici, pysb2jax
from amici.sim.jax import JAXProblem, ReturnValue, run_simulations
from amici.sim.sundials import (
    ExpData,
    SensitivityMethod,
    SensitivityOrder,
)
from amici.testing import TemporaryDirectoryWinSafe, skip_on_valgrind
from beartype import beartype
from numpy.testing import assert_allclose
from petab.v1.C import PREEQUILIBRATION_CONDITION_ID, SIMULATION_CONDITION_ID
from test_petab_objective import lotka_volterra  # noqa: F401

pysb = pytest.importorskip("pysb")

jax.config.update("jax_enable_x64", True)


ATOL_SIM = 1e-12
RTOL_SIM = 1e-12


@skip_on_valgrind
def test_code_printer_floats_roundtrip():
    """Generated code must recover the exact double.

    sympy's string printer emits only 15 significant digits, which perturbs
    the last bits and breaks exact comparisons in the generated model
    (SBML test suite case 00958: a parameter with value `pi` must test equal
    to `pi`).
    """
    printer = AmiciJaxCodePrinter()
    values = [sp.pi.evalf(), sp.Float(0.1), sp.Float(1) / 3, sp.Float(1e-17)]

    for value in values:
        assert float(printer.doprint(value)) == float(value)
        # `Max`/`Min` build their argument array themselves
        assert repr(float(value)) in printer.doprint(
            sp.Max(value, sp.Symbol("x"))
        )

    # parameter values are emitted without going through `doprint`
    assert _jnp_array_str(values) == "jnp.array([{}])".format(
        ", ".join(repr(float(value)) for value in values)
    )


@skip_on_valgrind
def test_conversion():
    pysb.SelfExporter.cleanup()  # reset pysb
    pysb.SelfExporter.do_export = True

    model = pysb.Model("conversion")
    a = pysb.Monomer("A", sites=["s"], site_states={"s": ["a", "b"]})
    pysb.Initial(a(s="a"), pysb.Parameter("aa0", 1.2))
    pysb.Rule("conv", a(s="a") >> a(s="b"), pysb.Parameter("kcat", 0.05))
    pysb.Observable("ab", a(s="b"))

    with TemporaryDirectoryWinSafe() as outdir:
        pysb2amici(model, outdir, verbose=True, observation_model=[MC("ab")])
        pysb2jax(model, outdir, verbose=True, observation_model=[MC("ab")])

        amici_module = import_model_module(
            module_name=model.name, module_path=outdir
        )
        jax_module = import_model_module(
            module_name=Path(outdir).stem, module_path=Path(outdir).parent
        )

        ts = tuple(np.linspace(0, 1, 10))
        p = jnp.stack((1.0, 0.1), axis=-1)
        k = tuple()
        _test_model(amici_module, jax_module, ts, p, k)


@skip_on_valgrind
@pytest.mark.filterwarnings(
    "ignore:Model does not contain any initial conditions"
)
def test_dimerization():
    pysb.SelfExporter.cleanup()  # reset pysb
    pysb.SelfExporter.do_export = True

    model = pysb.Model("dimerization")
    a = pysb.Monomer("A", sites=["b"])
    b = pysb.Monomer("B", sites=["a"])

    pysb.Rule(
        "turnover_a",
        a(b=None) | None,
        pysb.Parameter("kdeg_a", 10),
        pysb.Parameter("ksyn_a", 0.1),
    )
    pysb.Rule(
        "turnover_b",
        b(a=None) | None,
        pysb.Parameter("kdeg_b", 0.1),
        pysb.Parameter("ksyn_b", 10),
    )
    pysb.Rule(
        "dimer",
        a(b=None) + b(a=None) | a(b=1) % b(a=1),
        pysb.Parameter("kon", 1.0),
        pysb.Parameter("koff", 0.1),
    )

    pysb.Observable("a_obs", a())
    pysb.Observable("b_obs", b())

    with TemporaryDirectoryWinSafe() as outdir:
        pysb2amici(
            model,
            outdir,
            verbose=True,
            observation_model=[MC("a_obs"), MC("b_obs")],
            fixed_parameters=["ksyn_a", "ksyn_b"],
        )
        pysb2jax(
            model,
            outdir,
            observation_model=[MC("a_obs"), MC("b_obs")],
        )

        amici_module = import_model_module(
            module_name=model.name, module_path=outdir
        )
        jax_module = import_model_module(
            module_name=Path(outdir).stem, module_path=Path(outdir).parent
        )

        ts = tuple(np.linspace(0, 1, 10))
        p = jnp.stack((5, 0.5, 0.5, 0.5), axis=-1)
        k = (0.5, 5)
        _test_model(amici_module, jax_module, ts, p, k)


def _test_model(amici_module, jax_module, ts, p, k):
    amici_model = amici_module.get_model()

    amici_model.set_timepoints(np.asarray(ts, dtype=np.float64))
    sol_amici_ref = amici_model.simulate()

    jax_model = jax_module.Model()

    amici_model.set_free_parameters(np.asarray(p, dtype=np.float64))
    amici_model.set_fixed_parameters(np.asarray(k, dtype=np.float64))
    edata = ExpData(sol_amici_ref, 1.0, 1.0, 1)
    edata.free_parameters = amici_model.get_free_parameters()
    edata.fixed_parameters = amici_model.get_fixed_parameters()
    edata.pscale = amici_model.get_parameter_scale()
    amici_solver = amici_model.create_solver()
    amici_solver.set_sensitivity_method(SensitivityMethod.forward)
    amici_solver.set_sensitivity_order(SensitivityOrder.first)
    amici_solver.set_absolute_tolerance(ATOL_SIM)
    amici_solver.set_relative_tolerance(RTOL_SIM)
    rs_amici = amici_model.simulate(solver=amici_solver, edata=[edata])

    check_fields_jax(
        rs_amici,
        jax_model,
        amici_model.get_free_parameter_ids(),
        amici_model.get_fixed_parameter_ids(),
        edata,
        ["x", "y", "llh", "res", "x0"],
    )

    check_fields_jax(
        rs_amici,
        jax_model,
        amici_model.get_free_parameter_ids(),
        amici_model.get_fixed_parameter_ids(),
        edata,
        ["sllh", "sx0", "sx", "sres", "sy"],
        sensi_order=SensitivityOrder.first,
    )


def check_fields_jax(
    rs_amici,
    jax_model,
    free_parameter_ids,
    fixed_parameter_ids,
    edata,
    fields,
    sensi_order=SensitivityOrder.none,
):
    r_jax = dict()
    ts = np.array(edata.get_timepoints())
    my = np.array(edata.get_measurements()).reshape(len(ts), -1)
    ts = np.repeat(ts.reshape(-1, 1), my.shape[1], axis=1)
    iys = np.repeat(np.arange(my.shape[1]).reshape(1, -1), len(ts), axis=0)
    my = my.flatten()
    ts = ts.flatten()
    iys = iys.flatten()
    iy_trafos = np.zeros_like(iys)

    ts_dyn = ts
    ts_posteq = np.array([])

    par_dict = {
        **dict(zip(free_parameter_ids, edata.free_parameters)),
        **dict(zip(fixed_parameter_ids, edata.fixed_parameters)),
    }

    p = jnp.array([par_dict[par_id] for par_id in jax_model.parameter_ids])
    kwargs = {
        "ts_dyn": jnp.array(ts_dyn),
        "ts_posteq": jnp.array(ts_posteq),
        "my": jnp.array(my),
        "iys": jnp.array(iys),
        "ops": jnp.zeros((*my.shape[:2], 0)),
        "nps": jnp.zeros((*my.shape[:2], 0)),
        "iy_trafos": jnp.array(iy_trafos),
        "x_preeq": jnp.array([]),
        "solver": diffrax.Kvaerno5(),
        "controller": diffrax.PIDController(atol=1e-8, rtol=1e-8),
        "root_finder": optimistix.Newton(atol=ATOL_SIM, rtol=RTOL_SIM),
        "adjoint": diffrax.RecursiveCheckpointAdjoint(),
        "steady_state_event": diffrax.steady_state_event(),
        "max_steps": 2**8,  # max_steps
    }
    # Use beartype-wrapped unjitted version for type checking
    # (beartype cannot introspect jitted functions, so we wrap the unjitted version)
    fun = beartype(jax_model.simulate_condition_unjitted)

    for output in ["llh", "x0", "x", "y", "res"]:
        okwargs = kwargs | {
            "adjoint": diffrax.DirectAdjoint(),
            "max_steps": 2**8,
            "ret": ReturnValue[output],
        }
        if sensi_order == SensitivityOrder.none:
            r_jax[output] = fun(p, **okwargs)[0]
        if sensi_order == SensitivityOrder.first:
            if output == "llh":
                r_jax[f"s{output}"] = jax.grad(fun, has_aux=True)(p, **kwargs)[
                    0
                ]
            else:
                r_jax[f"s{output}"] = jax.jacfwd(fun, has_aux=True)(
                    p, **okwargs
                )[0]

    amici_par_idx = np.array(
        [
            jax_model.parameter_ids.index(par_id)
            for par_id in free_parameter_ids
        ]
    )

    for field in fields:
        for r_amici, r_jax in zip(rs_amici, [r_jax]):
            actual = r_jax[field]
            desired = r_amici[field]
            if field == "x":
                actual = actual[iys == 0, :]
            if field == "y":
                actual = np.stack(
                    [actual[iys == iy] for iy in sorted(np.unique(iys))],
                    axis=1,
                )
            elif field == "sllh":
                actual = actual[amici_par_idx]
            elif field == "sx":
                actual = actual[:, :, amici_par_idx]
                actual = np.permute_dims(actual[iys == 0, :, :], (0, 2, 1))
            elif field == "sy":
                actual = actual[:, amici_par_idx]
                actual = np.permute_dims(
                    np.stack(
                        [
                            actual[iys == iy, :]
                            for iy in sorted(np.unique(iys))
                        ],
                        axis=1,
                    ),
                    (0, 2, 1),
                )
            elif field == "sx0":
                actual = actual[:, amici_par_idx].T
            elif field == "sres":
                actual = actual[:, amici_par_idx]

            assert_allclose(
                actual=actual,
                desired=desired,
                atol=1e-5,
                rtol=1e-5,
                err_msg=f"field {field} does not match",
            )


def test_preequilibration_failure(lotka_volterra):  # noqa: F811
    petab_problem = lotka_volterra
    # oscillating system, preequilibation should fail when interaction is active

    with TemporaryDirectoryWinSafe(prefix="normal") as model_dir:
        jax_problem = import_petab_problem(
            petab_problem, jax=True, output_dir=model_dir
        )
        r = run_simulations(jax_problem)
        assert not np.isinf(r[0].item())
    petab_problem.measurement_df[PREEQUILIBRATION_CONDITION_ID] = (
        petab_problem.measurement_df[SIMULATION_CONDITION_ID]
    )
    with TemporaryDirectoryWinSafe(prefix="failure") as model_dir:
        jax_problem = import_petab_problem(
            petab_problem, jax=True, output_dir=model_dir
        )
        r = run_simulations(jax_problem)
        assert np.isinf(r[0].item())


@skip_on_valgrind
def test_serialisation(lotka_volterra):  # noqa: F811
    petab_problem = lotka_volterra
    with TemporaryDirectoryWinSafe(
        prefix=petab_problem.model.model_id
    ) as model_dir:
        jax_problem = import_petab_problem(
            petab_problem, jax=True, output_dir=model_dir
        )
        # change parameters to random values to test serialisation
        jax_problem.update_parameters(
            jax_problem.parameters
            + jr.normal(jr.PRNGKey(0), jax_problem.parameters.shape)
        )

        with TemporaryDirectoryWinSafe() as outdir:
            outdir = Path(outdir)
            jax_problem.save(outdir)
            jax_problem_loaded = JAXProblem.load(outdir)
            assert_allclose(
                jax_problem.parameters, jax_problem_loaded.parameters
            )


@skip_on_valgrind
def test_condition_table_initial_value_is_differentiable(tmp_path):
    """A parameter used as a species initial value via the condition table
    must stay a live function of ``JAXProblem.parameters``.

    Regression test for a bug where the reinitialisation value was resolved
    through the ``targets_map`` cached at construction time
    (see ``JAXProblem._state_reinitialisation_value``): the initial value was
    then frozen at the nominal parameter value, so ``update_parameters`` had
    no effect on it and its gradient silently leaked into the cache instead of
    ``grad.parameters``.
    """
    import equinox as eqx
    import petab.v1 as petab
    from petab.v1.models.sbml_model import SbmlModel

    problem = petab.Problem()
    problem.model = SbmlModel.from_antimony(
        "compartment_ = 1;\n"
        "species A in compartment_, B in compartment_;\n"
        "A = 3; B = 0;\n"
        "k1 = 0.8; k2 = 0.6;\n"
        "fwd: A -> B; k1 * A;\n"
        "rev: B -> A; k2 * B;\n"
    )
    # `a0` initialises species `A` via the condition table (not in the model)
    problem.add_parameter(
        "a0", estimate=True, nominal_value=2.0, scale="lin", lb=0.1, ub=10
    )
    problem.add_observable("obs_a", "A", noise_formula="0.5")
    problem.add_condition("c0", A="a0")
    problem.add_measurement("obs_a", "c0", 0.0, 0.7)
    problem.add_measurement("obs_a", "c0", 10.0, 0.1)

    jax_problem = import_petab_problem(
        problem, jax=True, output_dir=str(tmp_path)
    )
    ia = jax_problem.parameter_ids.index("a0")

    def llh(p):
        return run_simulations(jax_problem.update_parameters(p))[0]

    p0 = jax_problem.parameters
    # `a0` enters the likelihood only through A(0); updating it must change llh
    assert abs(float(llh(p0.at[ia].add(1.0))) - float(llh(p0))) > 1e-6, (
        "initial value is frozen w.r.t. update_parameters"
    )

    # autodiff w.r.t. `a0` (read from grad.parameters) must match finite diff
    eps = 1e-6
    fd = (float(llh(p0.at[ia].add(eps))) - float(llh(p0.at[ia].add(-eps)))) / (
        2 * eps
    )
    grad = eqx.filter_grad(lambda m: run_simulations(m)[0])(
        jax_problem.update_parameters(p0)
    )
    assert_allclose(float(grad.parameters[ia]), fd, rtol=1e-4, atol=1e-4)


@skip_on_valgrind
def test_condition_table_parameter_override_is_differentiable(tmp_path):
    """A model parameter mapped to an estimated parameter via the condition
    table (the standard PEtab pattern for condition-specific estimated
    parameters) must stay a live function of ``JAXProblem.parameters``.

    Regression test for the same construction-time freezing bug as
    ``test_condition_table_initial_value_is_differentiable``, on the parameter
    mapping path (``JAXProblem._map_experiment_model_parameter_value``).
    """
    import equinox as eqx
    import petab.v1 as petab
    from petab.v1.models.sbml_model import SbmlModel

    problem = petab.Problem()
    problem.model = SbmlModel.from_antimony(
        "compartment_ = 1;\n"
        "species A in compartment_, B in compartment_;\n"
        "A = 1; B = 0;\n"
        "k1 = 0.8; k2 = 0.6;\n"
        "fwd: A -> B; k1 * A;\n"
        "rev: B -> A; k2 * B;\n"
    )
    # the condition maps model parameter `k1` to the estimated `k1_c0`
    problem.add_parameter(
        "k1_c0", estimate=True, nominal_value=0.8, scale="lin", lb=0.1, ub=10
    )
    problem.add_observable("obs_b", "B", noise_formula="0.5")
    problem.add_condition("c0", k1="k1_c0")
    problem.add_measurement("obs_b", "c0", 1.0, 0.3)
    problem.add_measurement("obs_b", "c0", 5.0, 0.4)

    jax_problem = import_petab_problem(
        problem, jax=True, output_dir=str(tmp_path)
    )
    ik = jax_problem.parameter_ids.index("k1_c0")

    def llh(p):
        return run_simulations(jax_problem.update_parameters(p))[0]

    p0 = jax_problem.parameters
    assert abs(float(llh(p0.at[ik].add(0.3))) - float(llh(p0))) > 1e-6, (
        "condition-overridden parameter is frozen w.r.t. update_parameters"
    )

    eps = 1e-6
    fd = (float(llh(p0.at[ik].add(eps))) - float(llh(p0.at[ik].add(-eps)))) / (
        2 * eps
    )
    grad = eqx.filter_grad(lambda m: run_simulations(m)[0])(
        jax_problem.update_parameters(p0)
    )
    assert_allclose(float(grad.parameters[ik]), fd, rtol=1e-4, atol=1e-4)


@skip_on_valgrind
def test_petab_simulate_ragged_experiments(tmp_path):
    """``petab_simulate`` must handle experiments with different numbers of
    measurement timepoints.

    Regression test for ``_build_simulation_df_v2``: ``_get_measurements``
    pads every experiment's arrays to a common length, so an experiment with
    fewer timepoints has masked-out padding. If the padding mask is not
    applied consistently to the index and all columns, building the
    simulation DataFrame raises ``ValueError: arrays must all be same
    length`` (or leaks padded/duplicated indices).
    """
    import petab.v1 as petab
    from amici.sim.jax import petab_simulate
    from petab.v1.models.sbml_model import SbmlModel

    problem = petab.Problem()
    problem.model = SbmlModel.from_antimony(
        "compartment_ = 1;\n"
        "species A in compartment_, B in compartment_;\n"
        "A = 1; B = 0;\n"
        "k1 = 0.8; k2 = 0.6;\n"
        "fwd: A -> B; k1 * A;\n"
        "rev: B -> A; k2 * B;\n"
    )
    # `k1` is set per condition below, so only the free `k2` goes in the
    # parameter table
    problem.add_parameter(
        "k2", estimate=False, nominal_value=0.6, scale="lin", lb=0.1, ub=10
    )
    problem.add_observable("obs_b", "B", noise_formula="0.5")
    # two conditions -> two experiments, with DIFFERENT numbers of
    # timepoints so `_ts_masks` has genuine padding
    problem.add_condition("c0", k1=0.8)
    problem.add_condition("c1", k1=0.5)
    problem.add_measurement("obs_b", "c0", 1.0, 0.3)
    problem.add_measurement("obs_b", "c1", 1.0, 0.4)
    problem.add_measurement("obs_b", "c1", 5.0, 0.2)

    jax_problem = import_petab_problem(
        problem, jax=True, output_dir=str(tmp_path)
    )

    sim_df = petab_simulate(jax_problem)

    # exactly one simulated row per measurement (1 for c0, 2 for c1) -- no
    # length mismatch, no padded/duplicated rows, no missing simulations
    assert len(sim_df) == len(problem.measurement_df) == 3
    assert sim_df.index.is_unique
    assert sorted(sim_df[petab.TIME].tolist()) == [1.0, 1.0, 5.0]
    assert not sim_df[petab.SIMULATION].isna().any()


@skip_on_valgrind
def test_steady_state_event_no_recompile_across_conditions(
    tmp_path, monkeypatch
):
    """Simulating different conditions/parameters must not force JAX to
    recompile, as long as only numeric inputs (e.g. parameters) change.
    """
    import equinox as eqx
    from amici.importers.antimony import antimony2sbml
    from amici.importers.sbml import SbmlImporter
    from amici.sim.jax.model import JAXModel
    from amici.sim.jax.petab import (
        DEFAULT_CONTROLLER_SETTINGS,
        DEFAULT_ROOT_FINDER_SETTINGS,
        SteadyStateEvent,
    )

    ant_model = """
    model steady_state_recompile
        x' = -k * x
        x = 1
        k = 1
    end
    """
    sbml = antimony2sbml(ant_model)
    importer = SbmlImporter(sbml, from_file=False)
    importer.sbml2jax("steady_state_recompile", output_dir=tmp_path)
    module = amici._module_from_path(
        "steady_state_recompile", tmp_path / "__init__.py"
    )
    model = module.Model()

    def fresh_solver_kwargs():
        # construct brand new instances on every call, mimicking e.g. an
        # optimizer objective that rebuilds solver settings each iteration
        return dict(
            solver=diffrax.Kvaerno5(),
            controller=diffrax.PIDController(**DEFAULT_CONTROLLER_SETTINGS),
            root_finder=optimistix.Newton(**DEFAULT_ROOT_FINDER_SETTINGS),
            steady_state_event=SteadyStateEvent(),
        )

    conditions = [1.0, 2.5, 0.3]  # numeric-only differences between calls

    # Avoid order-dependent failures if prior tests already compiled these methods.
    jax.clear_caches()

    def patch_trace_counter(target, name):
        # eqx.debug.assert_max_traces/get_num_traces track calls to a plain
        # callable; since it isn't itself a descriptor, dispatch to it
        # manually so `self` is still bound correctly when called as
        # `self.<name>(...)`.
        wrapped = eqx.debug.assert_max_traces(
            getattr(target, name), max_traces=1
        )

        def dispatch(self, *args, **kwargs):
            return wrapped(self, *args, **kwargs)

        monkeypatch.setattr(target, name, dispatch)
        return wrapped

    ts = jnp.array([0.0, 1.0, 2.0])
    my = jnp.zeros_like(ts)
    iys = jnp.zeros_like(ts, dtype=int)
    iy_trafos = jnp.zeros_like(ts, dtype=int)

    simulate_traces = patch_trace_counter(
        JAXModel, "simulate_condition_unjitted"
    )
    for k_val in conditions:
        kwargs = fresh_solver_kwargs()
        model.simulate_condition(
            jnp.array([k_val]),
            ts,
            jnp.array([]),
            my,
            iys,
            iy_trafos,
            jnp.zeros((3, 0)),
            jnp.zeros((3, 0)),
            kwargs["solver"],
            kwargs["controller"],
            kwargs["root_finder"],
            diffrax.RecursiveCheckpointAdjoint(),
            kwargs["steady_state_event"],
            1000,
        )

    assert eqx.debug.get_num_traces(simulate_traces) == 1, (
        "simulate_condition was retraced across conditions with only "
        "numeric differences"
    )

    preeq_traces = patch_trace_counter(JAXModel, "_handle_t0_event")
    for k_val in conditions:
        kwargs = fresh_solver_kwargs()
        model.preequilibrate_condition(
            jnp.array([k_val]),
            jnp.array([]),
            jnp.array([]),
            jnp.array([]),
            kwargs["solver"],
            kwargs["controller"],
            kwargs["root_finder"],
            kwargs["steady_state_event"],
            1000,
        )

    assert eqx.debug.get_num_traces(preeq_traces) == 1, (
        "preequilibrate_condition was retraced across conditions with only "
        "numeric differences"
    )


@skip_on_valgrind
def test_time_dependent_discontinuity(tmp_path):
    """Models with time dependent discontinuities are handled."""

    from amici.importers.antimony import antimony2sbml
    from amici.importers.sbml import SbmlImporter
    from amici.sim.jax._simulation import solve
    from amici.sim.jax.petab import DEFAULT_CONTROLLER_SETTINGS

    ant_model = """
    model time_disc
        x' = piecewise(1, time - sin(time) - 1 < 0, 2)
        x = 0
    end
    """

    sbml = antimony2sbml(ant_model)
    importer = SbmlImporter(sbml, from_file=False)

    importer.sbml2jax("time_disc", output_dir=tmp_path)

    module = amici._module_from_path("time_disc", tmp_path / "__init__.py")
    model = module.Model()

    p = jnp.array([1.0])
    x0_full = model._x0(0.0, p)
    tcl = model._tcl(x0_full, p)
    x0 = model._x_solver(x0_full)
    ts = jnp.array([0.0, 1.0, 2.0])
    h = model._initialise_heaviside_variables(0.0, model._x_solver(x0), p, tcl)

    assert len(model._root_cond_fns()) > 0
    assert model._known_discs(p).size == 0

    ys, _, _ = solve(
        p,
        ts[0],
        ts,
        tcl,
        h,
        x0,
        jnp.ones_like(h),
        diffrax.Tsit5(),
        diffrax.PIDController(**DEFAULT_CONTROLLER_SETTINGS),
        optimistix.Newton(atol=1e-8, rtol=1e-8),
        1000,
        diffrax.DirectAdjoint(),
        diffrax.ODETerm(model._xdot),
        model._root_cond_fns(),
        model._root_cond_fn,
        model._delta_x,
        model._known_discs(p),
        model.observable_ids,
    )

    assert ys.shape[0] == ts.shape[0]


@skip_on_valgrind
def test_time_dependent_discontinuity_equilibration(tmp_path):
    """Time dependent discontinuities are handled during equilibration."""

    from amici.importers.antimony import antimony2sbml
    from amici.importers.sbml import SbmlImporter
    from amici.sim.jax._simulation import eq
    from amici.sim.jax.petab import DEFAULT_CONTROLLER_SETTINGS

    ant_model = """
    model time_disc_eq
        x' = piecewise(1, time - sin(time) - 1 < 0, -x)
        x = 0
    end
    """

    sbml = antimony2sbml(ant_model)
    importer = SbmlImporter(sbml, from_file=False)
    importer.sbml2jax("time_disc_eq", output_dir=tmp_path)

    module = amici._module_from_path("time_disc_eq", tmp_path / "__init__.py")
    model = module.Model()

    p = jnp.array([1.0])
    x0_full = model._x0(0.0, p)
    tcl = model._tcl(x0_full, p)
    x0 = model._x_solver(x0_full)
    h = model._initialise_heaviside_variables(0.0, model._x_solver(x0), p, tcl)

    assert len(model._root_cond_fns()) > 0
    assert model._known_discs(p).size == 0

    xs, _, _ = eq(
        p,
        tcl,
        h,
        x0,
        jnp.ones_like(h),
        diffrax.Tsit5(),
        diffrax.PIDController(**DEFAULT_CONTROLLER_SETTINGS),
        optimistix.Newton(atol=1e-8, rtol=1e-8),
        diffrax.steady_state_event(rtol=1e-8, atol=1e-8),
        diffrax.ODETerm(model._xdot),
        model._root_cond_fns(),
        model._root_cond_fn,
        model._delta_x,
        model._known_discs(p),
        1000,
    )

    assert_allclose(xs[0], 0.0, atol=1e-2)


@skip_on_valgrind
def test_explicit_discontinuity(tmp_path):
    """Explicit (time-triggered) discontinuities are emitted as tracked
    heaviside variables (not inlined ``jnp.select`` functions) and the solver
    is clipped onto the known discontinuity time.
    """
    from amici.importers.antimony import antimony2sbml
    from amici.importers.sbml import SbmlImporter
    from amici.sim.jax._simulation import solve
    from amici.sim.jax.petab import DEFAULT_CONTROLLER_SETTINGS

    # dx/dt = 1 for time > 5 else 0, x(0) = 0  =>  x(t) = max(t - 5, 0)
    ant_model = """
    model explicit_disc
        x' = piecewise(1, time > 5, 0)
        x = 0
    end
    """

    sbml = antimony2sbml(ant_model)
    importer = SbmlImporter(sbml, from_file=False)

    importer.sbml2jax("explicit_disc", output_dir=tmp_path)

    # the explicit heaviside must be a tracked variable, not an inlined
    # symbolic Heaviside function of time
    generated = (tmp_path / "__init__.py").read_text()
    assert "jnp.select" not in generated

    module = amici._module_from_path("explicit_disc", tmp_path / "__init__.py")
    model = module.Model()

    p = jnp.array([])
    x0_full = model._x0(0.0, p)
    tcl = model._tcl(x0_full, p)
    x0 = model._x_solver(x0_full)
    ts = jnp.array([0.0, 4.0, 5.0, 6.0, 10.0])
    h = model._initialise_heaviside_variables(0.0, x0, p, tcl)

    # explicit trigger time is a known discontinuity
    assert model._known_discs(p).size > 0

    ys, _, stats = solve(
        p,
        ts[0],
        ts,
        tcl,
        h,
        x0,
        jnp.ones_like(h),
        diffrax.Tsit5(),
        diffrax.PIDController(**DEFAULT_CONTROLLER_SETTINGS),
        optimistix.Newton(atol=1e-8, rtol=1e-8),
        1000,
        diffrax.DirectAdjoint(),
        diffrax.ODETerm(model._xdot),
        model._root_cond_fns(),
        model._root_cond_fn,
        model._delta_x,
        model._known_discs(p),
        model.observable_ids,
    )

    expected = jnp.clip(ts - 5.0, min=0.0)
    assert_allclose(ys.squeeze(), expected, atol=1e-6, rtol=1e-6)
    # clipping onto the known discontinuity avoids rejected steps
    assert stats["num_rejected_steps"] == 0


@skip_on_valgrind
def test_event_assignments_odd_root_count(tmp_path):
    """Event assignments are applied for every event, even when the total
    number of roots is odd.

    Regression test: ``_apply_event_assignments`` used to process roots in
    fixed windows of two (mirroring the root pairs generated for a
    persisted Heaviside variable), silently dropping the last root's
    assignment whenever a model's total root count was odd -- e.g. for any
    odd number of plain (non-Heaviside) events with assignments.
    """
    from amici.importers.antimony import antimony2sbml
    from amici.importers.sbml import SbmlImporter
    from amici.sim.jax._simulation import solve
    from amici.sim.jax.petab import DEFAULT_CONTROLLER_SETTINGS

    # three independent events, no Heaviside anywhere in the RHS -> three
    # roots total (odd), one per event
    ant_model = """
    model three_events
        x = 1
        y = 1
        z = 1
        x' = 0
        y' = 0
        z' = 0
        at time > 2: x = x + 10
        at time > 5: y = y + 100
        at time > 8: z = z + 1000
    end
    """

    sbml = antimony2sbml(ant_model)
    importer = SbmlImporter(sbml, from_file=False)

    importer.sbml2jax("three_events", output_dir=tmp_path)

    module = amici._module_from_path("three_events", tmp_path / "__init__.py")
    model = module.Model()

    p = jnp.array([])
    x0_full = model._x0(0.0, p)
    tcl = model._tcl(x0_full, p)
    x0 = model._x_solver(x0_full)
    ts = jnp.array([0.0, 1.0, 3.0, 6.0, 9.0, 10.0])
    h = model._initialise_heaviside_variables(0.0, x0, p, tcl)

    assert len(model._root_cond_fns()) == 3

    ys, _, _ = solve(
        p,
        ts[0],
        ts,
        tcl,
        h,
        x0,
        jnp.ones_like(h),
        diffrax.Tsit5(),
        diffrax.PIDController(**DEFAULT_CONTROLLER_SETTINGS),
        optimistix.Newton(atol=1e-8, rtol=1e-8),
        1000,
        diffrax.DirectAdjoint(),
        diffrax.ODETerm(model._xdot),
        model._root_cond_fns(),
        model._root_cond_fn,
        model._delta_x,
        model._known_discs(p),
        model.observable_ids,
    )

    expected = jnp.array(
        [
            [1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
            [11.0, 1.0, 1.0],
            [11.0, 101.0, 1.0],
            [11.0, 101.0, 1001.0],
            [11.0, 101.0, 1001.0],
        ]
    )
    assert_allclose(ys, expected, atol=1e-6, rtol=1e-6)
