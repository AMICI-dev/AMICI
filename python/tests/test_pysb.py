"""PYSB model tests"""
# flake8: noqa: F821

import importlib
import logging
import os

import pytest

pysb = pytest.importorskip("pysb")

from pathlib import Path

import amici
import numpy as np
import pysb.examples  # noqa: F811
import pytest
import sympy as sp
from amici import MeasurementChannel, import_model_module
from amici._symbolic.de_model_components import Event
from amici.importers.pysb import pysb2amici
from amici.importers.utils import amici_time_symbol, symbol_with_assumptions
from amici.sim.sundials import (
    ExpData,
    ParameterScaling,
    SensitivityOrder,
    SteadyStateSensitivityMode,
    parameter_scaling_from_int_vector,
    run_simulation,
)
from amici.sim.sundials.gradient_check import check_derivatives
from amici.testing import TemporaryDirectoryWinSafe, skip_on_valgrind
from numpy.testing import assert_allclose
from pysb.simulator import ScipyOdeSimulator


@skip_on_valgrind
def test_compare_to_sbml_import(
    pysb_example_presimulation_module, sbml_example_presimulation_module
):
    # -------------- PYSB -----------------

    model_pysb = pysb_example_presimulation_module.get_model()

    edata = get_data(model_pysb)

    rdata_pysb = get_results(model_pysb, edata)

    # -------------- SBML -----------------

    model_sbml = sbml_example_presimulation_module.get_model()

    rdata_sbml = get_results(model_sbml, edata)

    # check if preequilibration fixed parameters are correctly applied:
    for rdata, model, importer in zip(
        [rdata_sbml, rdata_pysb], [model_sbml, model_pysb], ["sbml", "pysb"]
    ):
        # check equilibrium fixed parameters
        assert_allclose(
            [sum(rdata["x_ss"][[1, 3]]), sum(rdata["x_ss"][[2, 4]])],
            edata.fixed_parameters_pre_equilibration,
            atol=1e-6,
            rtol=1e-6,
            err_msg=f"{importer} preequilibration",
        )

        # check equilibrium initial parameters
        assert_allclose(
            sum(rdata["x_ss"][[0, 3, 4, 5]]),
            model.get_free_parameter_by_name("PROT_0"),
            atol=1e-6,
            rtol=1e-6,
            err_msg=f"{importer} preequilibration",
        )

        # check reinitialization with fixed parameter after
        # presimulation
        assert_allclose(
            [rdata["x0"][1], rdata["x0"][2]],
            edata.fixed_parameters,
            atol=1e-6,
            rtol=1e-6,
            err_msg=f"{importer} presimulation",
        )

    skip_attrs = [
        "ptr",
        "preeq_t",
        "numsteps",
        "preeq_numsteps",
        "num_rhs_evals",
        "num_err_test_fails",
        "order",
        "J",
        "xdot",
        "preeq_wrms",
        "preeq_cpu_time",
        "cpu_time",
        "cpu_time_b",
        "cpu_time_total",
        "w",
    ]

    for field in rdata_pysb:
        if field in skip_attrs:
            continue

        if rdata_pysb[field] is None:
            assert rdata_sbml[field] is None, field
        elif rdata_sbml[field] is None:
            assert rdata_pysb[field] is None, field
        elif np.isnan(rdata_sbml[field]).all():
            assert np.isnan(rdata_pysb[field]).all(), field
        elif np.isnan(rdata_pysb[field]).all():
            assert np.isnan(rdata_sbml[field]).all(), field
        else:
            assert_allclose(
                rdata_sbml[field],
                rdata_pysb[field],
                atol=1e-6,
                rtol=1e-6,
                err_msg=field,
            )


pysb_models = [
    "tyson_oscillator",
    "robertson",
    "expression_observables",
    "bax_pore_sequential",
    "bax_pore",
    "bngwiki_egfr_simple",
    "bngwiki_enzymatic_cycle_mm",
    "bngwiki_simple",
    "earm_1_0",
    "earm_1_3",
    "move_connected",
    "michment",
    "kinase_cascade",
    "hello_pysb",
    "fricker_2010_apoptosis",
    "explicit",
    "fixed_initial",
    "localfunc",
]
custom_models = [
    "bngwiki_egfr_simple_deletemolecules",
]


@skip_on_valgrind
@pytest.mark.parametrize("example", pysb_models + custom_models)
def test_compare_to_pysb_simulation(example):
    atol = 1e-12
    rtol = 1e-10

    with amici.add_path(os.path.dirname(pysb.examples.__file__)):
        with amici.add_path(
            Path(__file__).parents[1] / "tests" / "pysb_test_models"
        ):
            # load example
            pysb.SelfExporter.cleanup()  # reset pysb
            pysb.SelfExporter.do_export = True

            module = importlib.import_module(example)
            pysb_model = module.model
            pysb_model.name = pysb_model.name.replace("pysb.examples.", "")
            # avoid naming clash for custom pysb models
            pysb_model.name += "_amici"

            # pysb part
            tspan = np.linspace(0, 100, 101)
            sim = ScipyOdeSimulator(
                pysb_model,
                tspan=tspan,
                integrator_options={"rtol": rtol, "atol": atol},
                compiler="python",
            )
            pysb_simres = sim.run()

            # amici part
            with TemporaryDirectoryWinSafe(prefix=pysb_model.name) as outdir:
                if pysb_model.name in ["move_connected_amici"]:
                    with pytest.raises(Exception):
                        pysb2amici(
                            pysb_model,
                            outdir,
                            verbose=logging.INFO,
                            compute_conservation_laws=True,
                        )
                    compute_conservation_laws = False
                else:
                    compute_conservation_laws = True

                pysb2amici(
                    pysb_model,
                    outdir,
                    verbose=logging.INFO,
                    compute_conservation_laws=compute_conservation_laws,
                    observation_model=list(
                        map(
                            MeasurementChannel,
                            pysb_model.observables.keys(),
                        )
                    ),
                )

                amici_model_module = import_model_module(
                    pysb_model.name, outdir
                )
                model_pysb = amici_model_module.get_model()
                model_pysb.set_timepoints(tspan)

                solver = model_pysb.create_solver()
                solver.set_max_steps(int(1e6))
                solver.set_absolute_tolerance(atol)
                solver.set_relative_tolerance(rtol)
                rdata = run_simulation(model_pysb, solver)

                # check agreement of species simulations
                assert_allclose(rdata["x"], pysb_simres.species, 1e-4, 1e-4)

                if example not in [
                    "fricker_2010_apoptosis",
                    "fixed_initial",
                    "bngwiki_egfr_simple_deletemolecules",
                ]:
                    if example in [
                        "tyson_oscillator",
                        "bax_pore_sequential",
                        "bax_pore",
                        "kinase_cascade",
                        "bngwiki_egfr_simple",
                        "bngwiki_enzymatic_cycle_mm",
                        "bngwiki_simple",
                    ]:
                        solver.set_absolute_tolerance(1e-14)
                        solver.set_relative_tolerance(1e-14)
                        epsilon = 1e-4
                    else:
                        solver.set_absolute_tolerance(1e-10)
                        solver.set_relative_tolerance(1e-10)
                        epsilon = 1e-3
                    model_pysb.set_parameter_scale(
                        parameter_scaling_from_int_vector(
                            [
                                ParameterScaling.log10
                                if p > 0
                                else ParameterScaling.none
                                for p in model_pysb.get_free_parameters()
                            ]
                        )
                    )
                    check_derivatives(
                        model_pysb,
                        solver=solver,
                        epsilon=epsilon,
                        rtol=1e-2,
                        atol=1e-2,
                        skip_zero_pars=True,
                    )


def get_data(model):
    solver = model.create_solver()
    model.set_timepoints(np.linspace(0, 60, 61))
    model.set_steady_state_sensitivity_mode(
        SteadyStateSensitivityMode.integrateIfNewtonFails
    )

    rdata = run_simulation(model, solver)
    edata = ExpData(rdata, 0.1, 0.0)
    edata.t_presim = 2
    edata.fixed_parameters = [10, 2]
    edata.fixed_parameters_presimulation = [3, 2]
    edata.fixed_parameters_pre_equilibration = [3, 0]
    edata.reinitialize_fixed_parameter_initial_states = True
    return edata


def get_results(model, edata):
    solver = model.create_solver()
    solver.set_sensitivity_order(SensitivityOrder.first)
    edata.reinitialize_fixed_parameter_initial_states = True
    model.set_timepoints(np.linspace(0, 60, 61))
    model.set_steady_state_sensitivity_mode(
        SteadyStateSensitivityMode.integrateIfNewtonFails
    )
    return run_simulation(model, solver, edata)


@skip_on_valgrind
def test_names_and_ids(pysb_example_presimulation_module):
    model_pysb = pysb_example_presimulation_module.get_model()
    expected = {
        "expression_ids": (
            "__s2",
            "__s1",
            "__s5",
            "pPROT",
            "tPROT",
            "initProt",
            "initDrug",
            "initKin",
            "pPROT_obs",
        ),
        "fixed_parameter_ids": ("DRUG_0", "KIN_0"),
        "fixed_parameter_names": ("DRUG_0", "KIN_0"),
        "observable_ids": ("pPROT_obs",),
        "observable_names": ("pPROT_obs",),
        "free_parameter_ids": (
            "PROT_0",
            "kon_prot_drug",
            "koff_prot_drug",
            "kon_prot_kin",
            "kphospho_prot_kin",
            "kdephospho_prot",
        ),
        "state_ids": ("__s0", "__s1", "__s2", "__s3", "__s4", "__s5"),
        "state_names": (
            "PROT(kin=None, drug=None, phospho='u')",
            "DRUG(bound=None)",
            "KIN(bound=None)",
            "DRUG(bound=1) % PROT(kin=None, drug=1, phospho='u')",
            "KIN(bound=1) % PROT(kin=1, drug=None, phospho='u')",
            "PROT(kin=None, drug=None, phospho='p')",
        ),
    }
    # Names and IDs are the same here
    expected["expression_names"] = expected["expression_ids"]
    expected["free_parameter_names"] = expected["free_parameter_ids"]

    for field_name, cur_expected in expected.items():
        actual = getattr(model_pysb, f"get_{field_name}")()
        assert actual == cur_expected


@skip_on_valgrind
def test_heaviside_and_special_symbols():
    pysb.SelfExporter.cleanup()  # reset pysb
    pysb.SelfExporter.do_export = True

    model = pysb.Model("piecewise_test")
    a = pysb.Monomer("A")
    pysb.Initial(a(), pysb.Parameter("a0"))
    pysb.Rule(
        "deg",
        a() >> None,
        pysb.Expression(
            "rate",
            sp.Piecewise((1, pysb.Observable("a", a()) < 1), (0.0, True)),
        ),
    )

    with TemporaryDirectoryWinSafe(prefix=model.name) as outdir:
        amici_model = pysb2amici(
            model,
            outdir,
            verbose=True,
            observation_model=[MeasurementChannel("a")],
        )

        assert amici_model.ne


@skip_on_valgrind
def test_energy():
    pysb.SelfExporter.cleanup()
    model_pysb = pysb.Model("energy")
    pysb.Monomer("A", ["a", "b"])
    pysb.Monomer("B", ["a"])
    pysb.Parameter("RT", 2)
    pysb.Parameter("A_0", 10)
    pysb.Parameter("AB_0", 10)
    pysb.Parameter("phi", 0.5)
    pysb.Expression("E_AAB_RT", -5 / RT)
    pysb.Expression("E0_AA_RT", -1 / RT)
    pysb.Rule(
        "A_dimerize",
        A(a=None) + A(a=None) | A(a=1) % A(a=1),
        phi,
        E0_AA_RT,
        energy=True,
    )
    pysb.EnergyPattern("epAAB", A(a=1) % A(a=1, b=2) % B(a=2), E_AAB_RT)
    pysb.Initial(A(a=None, b=None), A_0)
    pysb.Initial(A(a=None, b=1) % B(a=1), AB_0)

    with TemporaryDirectoryWinSafe(prefix=model_pysb.name) as outdir:
        pysb2amici(model_pysb, output_dir=outdir)

        model_module = import_model_module(
            module_name=model_pysb.name, module_path=outdir
        )
        amici_model = model_module.get_model()
        amici_model.set_timepoints(np.logspace(-4, 5, 10))
        solver = amici_model.create_solver()
        solver.set_relative_tolerance(1e-14)
        solver.set_absolute_tolerance(1e-14)

        check_derivatives(
            amici_model, solver=solver, epsilon=1e-4, rtol=1e-2, atol=1e-2
        )


@skip_on_valgrind
def test_pysb_event(tempdir):
    """Test adding events to PySB models."""
    pysb.SelfExporter.cleanup()  # reset pysb
    pysb.SelfExporter.do_export = True

    model = pysb.Model("pysb_event_test")
    a = pysb.Monomer("A")
    pysb.Initial(a(), pysb.Parameter("a0"))
    pysb.Rule("deg", a() >> None, pysb.Parameter("k", 1.0))

    events = [
        Event(
            symbol=symbol_with_assumptions("event1"),
            name="Event1",
            value=amici_time_symbol - 5,
            assignments={
                symbol_with_assumptions("__s0"): symbol_with_assumptions(
                    "__s0"
                )
                + 1000
            },
            use_values_from_trigger_time=False,
        )
    ]

    outdir = tempdir
    pysb2amici(
        model,
        outdir,
        verbose=True,
        observation_model=[MeasurementChannel("a")],
        compute_conservation_laws=False,
        _events=events,
    )

    model_module = import_model_module(
        module_name=model.name, module_path=outdir
    )
    amici_model = model_module.get_model()
    assert amici_model.ne
    amici_model.set_timepoints([0, 4, 5])

    np.testing.assert_allclose(
        amici_model.simulate().x, np.array([[0.0], [0.0], [1000.0]])
    )


@skip_on_valgrind
def test_pysb_reserved_names(tempdir):
    """Model quantities literally named one of AMICI's reserved
    array-parameter names (x, p, k, h, w, y) are renamed internally and
    have their original ids restored everywhere they're reported --
    mirrors SBML import's handling of the same collision class, see
    test_reserved_symbols.py."""
    pysb.SelfExporter.cleanup()  # reset pysb
    pysb.SelfExporter.do_export = True

    model = pysb.Model("pysb_reserved_names_test")
    a = pysb.Monomer("A")
    pysb.Initial(a(), pysb.Parameter("a0", 1.0))
    pysb.Rule("deg", a() >> None, pysb.Parameter("p", 0.5))
    # reference "p" here too, so its renamed symbol survives substitution
    # into another expression
    pysb.Expression("h", 2 * model.parameters["p"])
    pysb.Observable("y", a())

    outdir = tempdir
    pysb2amici(
        model,
        outdir,
        observation_model=[MeasurementChannel("y")],
        compute_conservation_laws=False,
    )

    model_module = import_model_module(
        module_name=model.name, module_path=outdir
    )
    amici_model = model_module.get_model()

    assert amici_model.get_observable_ids() == ("y",)
    assert "p" in amici_model.get_free_parameter_ids()
    assert "h" in amici_model.get_expression_ids()

    amici_model.set_timepoints([0, 1, 2])
    solver = amici_model.create_solver()
    rdata = run_simulation(amici_model, solver)
    assert rdata.status == amici.sim.sundials.AMICI_SUCCESS
    assert np.all(np.isfinite(rdata.x))


@skip_on_valgrind
def test_pysb_derived_name_collisions(tempdir):
    """Model quantities literally named after AMICI-internally-derived
    symbols (a conservation-law total `tcl_s{n}`, or an observable's
    measurement/sigma symbols `m{obs}`/`sigma_{obs}`) must not collide
    with them -- the derived symbol is disambiguated instead, with the
    user's own parameter keeping its original id. See
    https://github.com/AMICI-dev/AMICI/issues/3240 and the analogous SBML
    coverage in test_derived_symbol_collisions.py."""
    pysb.SelfExporter.cleanup()  # reset pysb
    pysb.SelfExporter.do_export = True

    model = pysb.Model("pysb_derived_name_collision_test")
    a = pysb.Monomer("A", ["b"])
    b = pysb.Monomer("B", ["a"])
    pysb.Parameter("A_0", 10)
    pysb.Parameter("B_0", 5)
    pysb.Parameter("kf", 1.0)
    pysb.Parameter("kr", 0.5)
    # decoys: whichever species index ends up being the conservation-law
    # target, one of these collides with its `tcl_s{n}` total-abundance
    # symbol
    for i in range(4):
        pysb.Parameter(f"tcl_s{i}", 42.0)
    # decoys colliding with the observable's derived measurement/sigma
    # symbols (`m{obs}` / `sigma_{obs}`)
    pysb.Parameter("mobs_a", 1.0)
    pysb.Parameter("sigma_obs_a", 1.0)

    pysb.Initial(a(b=None), model.parameters["A_0"])
    pysb.Initial(b(a=None), model.parameters["B_0"])
    pysb.Rule(
        "bind",
        a(b=None) + b(a=None) | a(b=1) % b(a=1),
        model.parameters["kf"],
        model.parameters["kr"],
    )
    pysb.Observable("obs_a", a())

    outdir = tempdir
    pysb2amici(
        model,
        outdir,
        observation_model=[MeasurementChannel("obs_a")],
        compute_conservation_laws=True,
    )

    model_module = import_model_module(
        module_name=model.name, module_path=outdir
    )
    amici_model = model_module.get_model()

    # the decoy parameters kept their own original ids
    free_ids = amici_model.get_free_parameter_ids()
    assert all(f"tcl_s{i}" in free_ids for i in range(4))
    assert "mobs_a" in free_ids
    assert "sigma_obs_a" in free_ids

    amici_model.set_timepoints([0, 1, 2, 5])
    solver = amici_model.create_solver()
    rdata = run_simulation(amici_model, solver)
    assert rdata.status == amici.sim.sundials.AMICI_SUCCESS
    assert np.all(np.isfinite(rdata.x))


def test_pysb_named_sigma():
    """Observables with a sigma given as the name of a `pysb.Expression`,
    mixed with observables using the default sigma.

    Sigmas must be ordered like the observables, and the sigma expression
    must remain a regular expression, since it may be used elsewhere in the
    model. See https://github.com/AMICI-dev/AMICI/issues/3538.
    """
    from amici.importers.pysb import ode_model_from_pysb_importer

    pysb.SelfExporter.cleanup()  # reset pysb
    pysb.SelfExporter.do_export = True

    model = pysb.Model("pysb_named_sigma_test")
    a = pysb.Monomer("A")
    b = pysb.Monomer("B")
    pysb.Initial(a(), pysb.Parameter("A_0", 1.0))
    pysb.Initial(b(), pysb.Parameter("B_0", 1.0))
    p_sd = pysb.Parameter("p_sd", 0.1)
    obs_a = pysb.Observable("obs_a", a())
    pysb.Observable("obs_b", b())
    sd_a = pysb.Expression("sd_a", 2 * p_sd * obs_a)
    pysb.Rule("deg_a", a() >> None, pysb.Parameter("k", 1.0))
    # the sigma expression is also used in the dynamics
    pysb.Rule("deg_b", b() >> None, sd_a)

    ode = ode_model_from_pysb_importer(
        model,
        observation_model=[
            MeasurementChannel("obs_a", sigma="sd_a"),
            MeasurementChannel("obs_b", noise_distribution="log-laplace"),
        ],
    )

    assert [str(y.get_sym()) for y in ode.observables()] == [
        "obs_a",
        "obs_b",
    ]
    assert [str(s.get_sym()) for s in ode.sigma_ys()] == [
        "sigma_obs_a",
        "sigma_obs_b",
    ]
    assert [s.get_val() for s in ode.sigma_ys()] == [
        symbol_with_assumptions("sd_a"),
        sp.Float(1.0),
    ]
    assert "sd_a" in {str(w.get_sym()) for w in ode.expressions()}
    assert ode.get_noise_distributions() == ["normal", "log-laplace"]
    assert list(ode.eq("sigmay")) == [
        2 * symbol_with_assumptions("p_sd") * symbol_with_assumptions("obs_a"),
        sp.Float(1.0),
    ]


def test_pysb_events_disable_conservation_laws(caplog):
    """Conservation laws are turned off for PySB models with events, and
    event assignments to non-states are rejected during import."""
    from amici.importers.pysb import ode_model_from_pysb_importer

    model = pysb.Model("pysb_event_cl_test", _export=False)
    model.add_component(pysb.Monomer("A", ["b"], _export=False))
    model.add_component(pysb.Monomer("B", ["a"], _export=False))
    model.add_component(pysb.Parameter("a0", 1.0, _export=False))
    model.add_component(pysb.Parameter("k", 1.0, _export=False))
    model.add_initial(
        pysb.Initial(
            model.monomers["A"](b=None), model.parameters["a0"], _export=False
        )
    )
    model.add_initial(
        pysb.Initial(
            model.monomers["B"](a=None), model.parameters["a0"], _export=False
        )
    )
    model.add_component(
        pysb.Rule(
            "binding",
            model.monomers["A"](b=None) + model.monomers["B"](a=None)
            >> model.monomers["A"](b=1) % model.monomers["B"](a=1),
            model.parameters["k"],
            _export=False,
        )
    )

    def event(target):
        return Event(
            symbol=symbol_with_assumptions("event1"),
            name="event1",
            value=amici_time_symbol - 5,
            assignments={target: 2},
            use_values_from_trigger_time=False,
        )

    # without events, total A and total B are conserved
    ode = ode_model_from_pysb_importer(model)
    assert ode.num_cons_law() == 2

    with caplog.at_level(logging.WARNING, logger="amici.importers.pysb"):
        ode = ode_model_from_pysb_importer(
            model, _events=[event(symbol_with_assumptions("__s0"))]
        )
    assert ode.num_cons_law() == 0
    assert "Conservation laws are currently not supported" in caplog.text

    with pytest.raises(ValueError, match=r"\['foo'\].*not \(solver\) state"):
        ode_model_from_pysb_importer(
            model, _events=[event(symbol_with_assumptions("foo"))]
        )


@pytest.mark.parametrize("first", ["stau", "drootdp", "dtaudp", "drootdtcl"])
def test_conservation_law_derivatives_independent_of_order(first):
    """Derivatives involving conservation laws don't depend on the order in
    which they are computed."""
    from amici.importers.pysb import ode_model_from_pysb_importer

    def build():
        model = pysb.Model("lock_order", _export=False)
        r = pysb.Monomer("R", ["s"], {"s": ["u", "p"]}, _export=False)
        c = pysb.Monomer("C", _export=False)
        r_0 = pysb.Parameter("R_0", 2.0, _export=False)
        k1 = pysb.Parameter("k1", 1.0, _export=False)
        rp = pysb.Observable("Rp", r(s="p"), _export=False)
        # root function depending on `w`
        rate = pysb.Expression(
            "rate", sp.Piecewise((1, rp > 1.5), (0, True)), _export=False
        )
        for component in (r, c, r_0, k1, rp, rate):
            model.add_component(component)
        model.add_initial(pysb.Initial(r(s="u"), r_0, _export=False))
        for rule in (
            pysb.Rule("conv", r(s="u") >> r(s="p"), k1, _export=False),
            pysb.Rule("prod", None >> c(), rate, _export=False),
        ):
            model.add_component(rule)
        de_model = ode_model_from_pysb_importer(
            model, compute_conservation_laws=True
        )
        assert de_model.num_cons_law() == 1
        return de_model

    names = ("dwdp", "drootdp", "sroot", "stau", "dtaudp", "drootdtcl")
    ref = build()
    expected = {name: ref.eq(name) for name in names}

    # The conserved total of `R` enters the root functions via `Rp`, so
    #  `drootdp` contains the total's sensitivities, each exactly once.
    stcl = {
        str(s): s
        for s in expected["drootdp"].free_symbols
        if str(s).startswith("stcl_")
    }
    assert set(stcl) == {"stcl_s1__R_0", "stcl_s1__k1"}
    stcl_r0, stcl_k1 = stcl["stcl_s1__R_0"], stcl["stcl_s1__k1"]
    assert expected["drootdp"] == sp.Matrix(
        [[-stcl_r0, -stcl_k1], [stcl_r0, stcl_k1]]
    )
    for stau in expected["stau"]:
        assert set(stcl.values()) <= stau.free_symbols

    de_model = build()
    de_model.eq(first)
    for name in names:
        actual = de_model.eq(name)
        if isinstance(actual, sp.MatrixBase):
            actual, ref_val = [actual], [expected[name]]
        else:
            # per-event lists for `stau` and `dtaudp`
            ref_val = expected[name]
        for a, b in zip(actual, ref_val, strict=True):
            assert (a - b).is_zero_matrix, name
