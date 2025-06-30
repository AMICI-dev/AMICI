import pytest
import amici
from pathlib import Path

pytest.importorskip("jax")
import amici.jax

import jax.numpy as jnp
import jax.random as jr
import jax
import diffrax
import optimistix
import numpy as np
from beartype import beartype
from petab.v1.C import PREEQUILIBRATION_CONDITION_ID, SIMULATION_CONDITION_ID

from amici.pysb_import import pysb2amici, pysb2jax
from amici.testing import TemporaryDirectoryWinSafe, skip_on_valgrind
from amici.petab.petab_import import import_petab_problem
from amici.jax import JAXProblem, ReturnValue, run_simulations
from numpy.testing import assert_allclose
from test_petab_objective import lotka_volterra  # noqa: F401

pysb = pytest.importorskip("pysb")

jax.config.update("jax_enable_x64", True)


ATOL_SIM = 1e-12
RTOL_SIM = 1e-12


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
        pysb2amici(model, outdir, verbose=True, observables=["ab"])
        pysb2jax(model, outdir, verbose=True, observables=["ab"])

        amici_module = amici.import_model_module(
            module_name=model.name, module_path=outdir
        )
        jax_module = amici.import_model_module(
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
            observables=["a_obs", "b_obs"],
            constant_parameters=["ksyn_a", "ksyn_b"],
        )
        pysb2jax(
            model,
            outdir,
            observables=["a_obs", "b_obs"],
        )

        amici_module = amici.import_model_module(
            module_name=model.name, module_path=outdir
        )
        jax_module = amici.import_model_module(
            module_name=Path(outdir).stem, module_path=Path(outdir).parent
        )

        ts = tuple(np.linspace(0, 1, 10))
        p = jnp.stack((5, 0.5, 0.5, 0.5), axis=-1)
        k = (0.5, 5)
        _test_model(amici_module, jax_module, ts, p, k)


def _test_model(amici_module, jax_module, ts, p, k):
    amici_model = amici_module.getModel()

    amici_model.setTimepoints(np.asarray(ts, dtype=np.float64))
    sol_amici_ref = amici.runAmiciSimulation(
        amici_model, amici_model.getSolver()
    )

    jax_model = jax_module.Model()

    amici_model.setParameters(np.asarray(p, dtype=np.float64))
    amici_model.setFixedParameters(np.asarray(k, dtype=np.float64))
    edata = amici.ExpData(sol_amici_ref, 1.0, 1.0)
    edata.parameters = amici_model.getParameters()
    edata.fixedParameters = amici_model.getFixedParameters()
    edata.pscale = amici_model.getParameterScale()
    amici_solver = amici_model.getSolver()
    amici_solver.setSensitivityMethod(amici.SensitivityMethod.forward)
    amici_solver.setSensitivityOrder(amici.SensitivityOrder.first)
    amici_solver.setAbsoluteTolerance(ATOL_SIM)
    amici_solver.setRelativeTolerance(RTOL_SIM)
    rs_amici = amici.runAmiciSimulations(amici_model, amici_solver, [edata])

    check_fields_jax(
        rs_amici,
        jax_model,
        amici_model.getParameterIds(),
        amici_model.getFixedParameterIds(),
        edata,
        ["x", "y", "llh", "res", "x0"],
    )

    check_fields_jax(
        rs_amici,
        jax_model,
        amici_model.getParameterIds(),
        amici_model.getFixedParameterIds(),
        edata,
        ["sllh", "sx0", "sx", "sres", "sy"],
        sensi_order=amici.SensitivityOrder.first,
    )


def check_fields_jax(
    rs_amici,
    jax_model,
    parameter_ids,
    fixed_parameter_ids,
    edata,
    fields,
    sensi_order=amici.SensitivityOrder.none,
):
    r_jax = dict()
    ts = np.array(edata.getTimepoints())
    my = np.array(edata.getObservedData()).reshape(len(ts), -1)
    ts = np.repeat(ts.reshape(-1, 1), my.shape[1], axis=1)
    iys = np.repeat(np.arange(my.shape[1]).reshape(1, -1), len(ts), axis=0)
    my = my.flatten()
    ts = ts.flatten()
    iys = iys.flatten()
    iy_trafos = np.zeros_like(iys)

    ts_dyn = ts
    ts_posteq = np.array([])

    par_dict = {
        **dict(zip(parameter_ids, edata.parameters)),
        **dict(zip(fixed_parameter_ids, edata.fixedParameters)),
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
        "controller": diffrax.PIDController(atol=ATOL_SIM, rtol=RTOL_SIM),
        "root_finder": optimistix.Newton(atol=ATOL_SIM, rtol=RTOL_SIM),
        "adjoint": diffrax.RecursiveCheckpointAdjoint(),
        "steady_state_event": diffrax.steady_state_event(),
        "max_steps": 2**8,  # max_steps
    }
    fun = beartype(jax_model.simulate_condition)

    for output in ["llh", "x0", "x", "y", "res"]:
        okwargs = kwargs | {
            "adjoint": diffrax.DirectAdjoint(),
            "max_steps": 2**8,
            "ret": ReturnValue[output],
        }
        if sensi_order == amici.SensitivityOrder.none:
            r_jax[output] = fun(p, **okwargs)[0]
        if sensi_order == amici.SensitivityOrder.first:
            if output == "llh":
                r_jax[f"s{output}"] = jax.grad(fun, has_aux=True)(p, **kwargs)[
                    0
                ]
            else:
                r_jax[f"s{output}"] = jax.jacfwd(fun, has_aux=True)(
                    p, **okwargs
                )[0]

    amici_par_idx = np.array(
        [jax_model.parameter_ids.index(par_id) for par_id in parameter_ids]
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
        jax_model = import_petab_problem(
            petab_problem, jax=True, model_output_dir=model_dir
        )
        jax_problem = JAXProblem(jax_model, petab_problem)
        r = run_simulations(jax_problem)
        assert not np.isinf(r[0].item())
    petab_problem.measurement_df[PREEQUILIBRATION_CONDITION_ID] = (
        petab_problem.measurement_df[SIMULATION_CONDITION_ID]
    )
    with TemporaryDirectoryWinSafe(prefix="failure") as model_dir:
        jax_model = import_petab_problem(
            petab_problem, jax=True, model_output_dir=model_dir
        )
        jax_problem = JAXProblem(jax_model, petab_problem)
        r = run_simulations(jax_problem)
        assert np.isinf(r[0].item())


@skip_on_valgrind
def test_serialisation(lotka_volterra):  # noqa: F811
    petab_problem = lotka_volterra
    with TemporaryDirectoryWinSafe(
        prefix=petab_problem.model.model_id
    ) as model_dir:
        jax_model = import_petab_problem(
            petab_problem, jax=True, model_output_dir=model_dir
        )
        jax_problem = JAXProblem(jax_model, petab_problem)
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
def test_time_dependent_discontinuity(tmp_path):
    """Models with time dependent discontinuities are handled."""

    from amici.antimony_import import antimony2sbml
    from amici.sbml_import import SbmlImporter
    from amici.jax.petab import DEFAULT_CONTROLLER_SETTINGS
    from amici.jax._simulation import solve

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
        ts,
        tcl,
        h,
        x0,
        diffrax.Tsit5(),
        diffrax.PIDController(**DEFAULT_CONTROLLER_SETTINGS),
        optimistix.Newton(atol=1e-8, rtol=1e-8),
        1000,
        diffrax.DirectAdjoint(),
        diffrax.ODETerm(model._xdot),
        model._root_cond_fns(),
        model._root_cond_fn,
        model._known_discs(p),
    )

    assert ys.shape[0] == ts.shape[0]


@skip_on_valgrind
def test_time_dependent_discontinuity_equilibration(tmp_path):
    """Time dependent discontinuities are handled during equilibration."""

    from amici.antimony_import import antimony2sbml
    from amici.sbml_import import SbmlImporter
    from amici.jax.petab import DEFAULT_CONTROLLER_SETTINGS
    from amici.jax._simulation import eq

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
        diffrax.Tsit5(),
        diffrax.PIDController(**DEFAULT_CONTROLLER_SETTINGS),
        optimistix.Newton(atol=1e-8, rtol=1e-8),
        diffrax.steady_state_event(rtol=1e-8, atol=1e-8),
        diffrax.ODETerm(model._xdot),
        model._root_cond_fns(),
        model._root_cond_fn,
        model._known_discs(p),
        1000,
    )

    assert_allclose(xs[0], 0.0, atol=1e-2)
