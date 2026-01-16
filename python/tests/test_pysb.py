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
    pysb.Rule("deg", a() >> None, pysb.Parameter("kk", 1.0))

    events = [
        Event(
            # note that unlike for SBML import, we must omit the real=True here
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
