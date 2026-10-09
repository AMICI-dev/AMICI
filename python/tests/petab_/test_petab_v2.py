import copy

import pytest
from amici.importers.petab import (
    PetabImporter,
    flatten_timepoint_specific_output_overrides,
    has_timepoint_specific_overrides,
    unflatten_simulation_df,
)
from amici.sim.sundials import SensitivityOrder
from amici.testing import TemporaryDirectoryWinSafe
from petab.v2 import C, Problem
from petab.v2.models.sbml_model import SbmlModel


def test_problem_has_timepoint_specific_overrides():
    """Test Problem.measurement_table_has_timepoint_specific_mappings"""
    problem = Problem()
    problem.add_measurement(
        obs_id="obs1",
        time=1.0,
        measurement=0.1,
        observable_parameters=["obsParOverride"],
    )
    problem.add_measurement(
        obs_id="obs1",
        time=1.0,
        measurement=0.2,
        observable_parameters=["obsParOverride2"],
    )
    assert has_timepoint_specific_overrides(problem) is True

    # different observables
    problem.measurement_tables[0].measurements[1].observable_id = "obs2"
    assert has_timepoint_specific_overrides(problem) is False

    # mixed numeric string
    problem.measurement_tables[0].measurements[1].observable_id = "obs1"
    problem.measurement_tables[0].measurements[1].observable_parameters = [
        "obsParOverride"
    ]
    assert has_timepoint_specific_overrides(problem) is False

    # different numeric values
    problem.measurement_tables[0].measurements[1].noise_parameters = [2.0]
    assert has_timepoint_specific_overrides(problem) is True
    assert (
        has_timepoint_specific_overrides(
            problem, ignore_scalar_numeric_noise_parameters=True
        )
        is False
    )


def test_flatten_timepoint_specific_output_overrides():
    """Test flatten_timepoint_specific_output_overrides"""
    problem = Problem()
    problem.model = SbmlModel.from_antimony("""x = 1""")
    for par_id in (
        "noiseParOverride1",
        "noiseParOverride2",
        "obsParOverride1",
        "obsParOverride2",
    ):
        problem.add_parameter(par_id, estimate=False, nominal_value=1)

    problem_expected = copy.deepcopy(problem)

    problem.add_observable(
        "obs1",
        formula="observableParameter1_obs1 + observableParameter2_obs1",
        noise_formula=(
            "(observableParameter1_obs1 + "
            "observableParameter2_obs1) * noiseParameter1_obs1"
        ),
        observable_placeholders=[
            "observableParameter1_obs1",
            "observableParameter2_obs1",
        ],
        noise_placeholders=["noiseParameter1_obs1"],
    )
    problem.add_observable("obs2", formula="x", noise_formula="1")

    # new observable IDs
    #  (obs${i_obs}_${i_obsParOverride}_${i_noiseParOverride})
    obs1_1_1 = "obs1__obsParOverride1_1_00000000000000__noiseParOverride1"
    obs1_2_1 = "obs1__obsParOverride2_1_00000000000000__noiseParOverride1"
    obs1_2_2 = "obs1__obsParOverride2_1_00000000000000__noiseParOverride2"

    for obs_id in (obs1_1_1, obs1_2_1, obs1_2_2):
        problem_expected.add_observable(
            obs_id,
            formula=(
                f"observableParameter1_{obs_id} "
                f"+ observableParameter2_{obs_id}"
            ),
            noise_formula=(
                f"(observableParameter1_{obs_id} + "
                f"observableParameter2_{obs_id}) "
                f"* noiseParameter1_{obs_id}"
            ),
            observable_placeholders=[
                f"observableParameter1_{obs_id}",
                f"observableParameter2_{obs_id}",
            ],
            noise_placeholders=[f"noiseParameter1_{obs_id}"],
        )

    problem_expected.add_observable(
        "obs2",
        formula="x",
        noise_formula="1",
    )

    # Measurement table with timepoint-specific overrides
    problem.add_measurement(
        obs_id="obs1",
        time=1.0,
        measurement=0.1,
        observable_parameters=["obsParOverride1", "1.0"],
        noise_parameters=["noiseParOverride1"],
    )
    problem.add_measurement(
        obs_id="obs1",
        time=1.0,
        measurement=0.1,
        observable_parameters=["obsParOverride2", "1.0"],
        noise_parameters=["noiseParOverride1"],
    )
    problem.add_measurement(
        obs_id="obs1",
        time=2.0,
        measurement=0.1,
        observable_parameters=["obsParOverride2", "1.0"],
        noise_parameters=["noiseParOverride2"],
    )
    problem.add_measurement(
        obs_id="obs1",
        time=2.0,
        measurement=0.1,
        observable_parameters=["obsParOverride2", "1.0"],
        noise_parameters=["noiseParOverride2"],
    )
    problem.add_measurement(obs_id="obs2", time=3.0, measurement=0.1)

    problem_expected.add_measurement(
        obs_id=obs1_1_1,
        time=1.0,
        measurement=0.1,
        observable_parameters=["obsParOverride1", "1.0"],
        noise_parameters=["noiseParOverride1"],
    )
    problem_expected.add_measurement(
        obs_id=obs1_2_1,
        time=1.0,
        measurement=0.1,
        observable_parameters=["obsParOverride2", "1.0"],
        noise_parameters=["noiseParOverride1"],
    )
    problem_expected.add_measurement(
        obs_id=obs1_2_2,
        time=2.0,
        measurement=0.1,
        observable_parameters=["obsParOverride2", "1.0"],
        noise_parameters=["noiseParOverride2"],
    )
    problem_expected.add_measurement(
        obs_id=obs1_2_2,
        time=2.0,
        measurement=0.1,
        observable_parameters=["obsParOverride2", "1.0"],
        noise_parameters=["noiseParOverride2"],
    )
    problem_expected.add_measurement(obs_id="obs2", time=3.0, measurement=0.1)

    problem.assert_valid()
    unflattened_problem = copy.deepcopy(problem)
    problem_expected.assert_valid()

    # Ensure having timepoint-specific overrides
    assert has_timepoint_specific_overrides(problem) is True
    assert has_timepoint_specific_overrides(problem_expected) is False

    flatten_timepoint_specific_output_overrides(problem)

    # Timepoint-specific overrides should be gone now
    assert has_timepoint_specific_overrides(problem) is False

    assert problem_expected.observables == problem.observables
    assert problem_expected.measurements == problem.measurements
    problem.assert_valid()

    simulation_df = copy.deepcopy(problem.measurement_df)
    simulation_df.rename(columns={C.MEASUREMENT: C.SIMULATION})
    unflattened_simulation_df = unflatten_simulation_df(
        simulation_df=simulation_df,
        petab_problem=unflattened_problem,
    )
    # The unflattened simulation dataframe has the original observable IDs.
    assert (
        unflattened_simulation_df[C.OBSERVABLE_ID] == ["obs1"] * 4 + ["obs2"]
    ).all()


def test_flatten_timepoint_specific_output_overrides_special_cases():
    """Test flatten_timepoint_specific_output_overrides
    for special cases:
    * no observable parameters
    """
    problem = Problem()
    problem.model = SbmlModel.from_antimony("""species1 = 1""")
    for p in ("noiseParOverride2", "noiseParOverride1"):
        problem.add_parameter(p, estimate=False, nominal_value=1)
    problem_expected = copy.deepcopy(problem)
    problem.add_observable(
        "obs1",
        formula="species1",
        noise_formula="noiseParameter1_obs1",
        noise_placeholders=["noiseParameter1_obs1"],
    )

    problem_expected.add_observable(
        "obs1__noiseParOverride1",
        formula="species1",
        noise_formula="noiseParameter1_obs1__noiseParOverride1",
        noise_placeholders=["noiseParameter1_obs1__noiseParOverride1"],
    )
    problem_expected.add_observable(
        "obs1__noiseParOverride2",
        formula="species1",
        noise_formula="noiseParameter1_obs1__noiseParOverride2",
        noise_placeholders=["noiseParameter1_obs1__noiseParOverride2"],
    )

    # Measurement table with timepoint-specific overrides
    problem.add_measurement(
        "obs1",
        time=1.0,
        measurement=0.1,
        noise_parameters=["noiseParOverride1"],
    )
    problem.add_measurement(
        "obs1",
        time=1.0,
        measurement=0.1,
        noise_parameters=["noiseParOverride1"],
    )
    problem.add_measurement(
        "obs1",
        time=2.0,
        measurement=0.1,
        noise_parameters=["noiseParOverride2"],
    )
    problem.add_measurement(
        "obs1",
        time=2.0,
        measurement=0.1,
        noise_parameters=["noiseParOverride2"],
    )

    problem_expected.add_measurement(
        "obs1__noiseParOverride1",
        time=1.0,
        measurement=0.1,
        noise_parameters=["noiseParOverride1"],
    )
    problem_expected.add_measurement(
        "obs1__noiseParOverride1",
        time=1.0,
        measurement=0.1,
        noise_parameters=["noiseParOverride1"],
    )
    problem_expected.add_measurement(
        "obs1__noiseParOverride2",
        time=2.0,
        measurement=0.1,
        noise_parameters=["noiseParOverride2"],
    )
    problem_expected.add_measurement(
        "obs1__noiseParOverride2",
        time=2.0,
        measurement=0.1,
        noise_parameters=["noiseParOverride2"],
    )

    problem.assert_valid()
    problem_expected.assert_valid()

    # Ensure having timepoint-specific overrides
    assert has_timepoint_specific_overrides(problem) is True

    flatten_timepoint_specific_output_overrides(problem)

    # Timepoint-specific overrides should be gone now
    assert has_timepoint_specific_overrides(problem) is False

    assert problem_expected.observables == problem.observables
    assert problem_expected.measurements == problem.measurements
    problem.assert_valid()


def test_petab_simulator_deepcopy_and_pickle():
    """Test that PetabImporter can be deep-copied"""
    problem = Problem()
    problem.model = SbmlModel.from_antimony("x = 1; x' = k;")
    problem.add_parameter("k", nominal_value=1.0, estimate=True, lb=0, ub=10)
    problem.add_observable("obs1", "x", noise_formula="1")
    for i in range(5):
        problem.add_measurement("obs1", time=i, measurement=2 * i)

    pi = PetabImporter(problem)
    ps = pi.create_simulator(force_import=False)
    ps.solver.set_sensitivity_order(SensitivityOrder.none)

    ps_copy = copy.deepcopy(ps)

    assert ps.simulate({"k": 2}).llh == ps_copy.simulate({"k": 2}).llh

    ps.solver.set_sensitivity_order(SensitivityOrder.first)
    assert (
        ps.solver.get_sensitivity_order()
        != ps_copy.solver.get_sensitivity_order()
    )

    import pickle

    ps_pickle = pickle.loads(pickle.dumps(ps))
    assert ps.simulate({"k": 2}).llh == ps_pickle.simulate({"k": 2}).llh


def _residual_test_problem() -> Problem:
    """A small two-experiment problem for residual (sensitivity) tests.

    All three observable placeholders are mapped to the same estimated
    parameter (``scale``), so that several model parameters map to a single
    problem parameter -- for ``obs_a`` even within a single residual. The
    noise formulas are parameter-independent, such that the least-squares
    identities ``FIM == sres.T @ sres`` and ``sllh == -res @ sres`` hold.
    """
    problem = Problem()
    problem.model = SbmlModel.from_antimony(
        "xa = 1; xb = 2; xa' = -k_decay * xa; xb' = -0.3 * xb; k_decay = 0.5;"
    )
    problem.add_parameter(
        "k_decay", nominal_value=0.5, estimate=True, lb=0.01, ub=10
    )
    problem.add_parameter(
        "scale", nominal_value=2.0, estimate=True, lb=0.01, ub=10
    )
    problem.add_observable(
        "obs_a",
        formula=(
            "observableParameter1_obs_a * xa + observableParameter2_obs_a"
        ),
        noise_formula="0.5",
        observable_placeholders=[
            "observableParameter1_obs_a",
            "observableParameter2_obs_a",
        ],
    )
    problem.add_observable(
        "obs_b",
        formula="observableParameter1_obs_b * xb",
        noise_formula="0.7",
        observable_placeholders=["observableParameter1_obs_b"],
    )
    # two experiments with different initial conditions, plus one experiment
    #  without any measurements, which must not contribute any residuals
    problem.add_condition("cond_lo", xa=1.0)
    problem.add_condition("cond_hi", xa=4.0)
    problem.add_experiment("exp_lo", 0, "cond_lo")
    problem.add_experiment("exp_hi", 0, "cond_hi")
    problem.add_experiment("exp_no_measurements", 0, "cond_lo")

    for experiment_id in ("exp_lo", "exp_hi"):
        for t in (0.0, 1.0, 2.0):
            problem.add_measurement(
                "obs_a",
                experiment_id=experiment_id,
                time=t,
                measurement=1.9 - 0.1 * t,
                observable_parameters=["scale", "scale"],
            )
            problem.add_measurement(
                "obs_b",
                experiment_id=experiment_id,
                time=t,
                measurement=3.8 - 0.2 * t,
                observable_parameters="scale",
            )
    problem.assert_valid()
    return problem


@pytest.fixture(scope="module")
def residual_test_importer() -> PetabImporter:
    """Importer for the model of :func:`_residual_test_problem`.

    The model is generated once for all tests using this fixture, into a
    temporary directory of its own -- so the tests neither depend on nor
    interfere with any other model module, also when run concurrently or
    distributed across processes. Each test creates its own simulator (and
    thus its own model and solver instance) from the importer.
    """
    with TemporaryDirectoryWinSafe(prefix="petab_v2_residuals_") as output_dir:
        pi = PetabImporter(
            _residual_test_problem(),
            module_name="test_petab_v2_residuals",
            output_dir=output_dir,
            verbose=False,
        )
        pi.import_module(force_import=True)
        yield pi


def test_aggregated_residual_sensitivities(residual_test_importer):
    """Aggregated residual sensitivities are consistent with the aggregated
    log-likelihood sensitivities, the FIM, and finite differences."""
    import numpy as np
    from amici.sim.sundials import SensitivityMethod

    problem = _residual_test_problem()
    ps = residual_test_importer.create_simulator()
    ps.solver.set_absolute_tolerance(1e-14)
    ps.solver.set_relative_tolerance(1e-12)

    x_ids = problem.x_free_ids
    x = {"k_decay": 0.4, "scale": 1.5}

    ps.solver.set_sensitivity_method(SensitivityMethod.forward)
    ps.solver.set_sensitivity_order(SensitivityOrder.first)
    result = ps.simulate(x)

    res = result.res
    sres = result.sres
    # 2 experiments (the third one has no measurements)
    #  x 3 timepoints x 2 observables
    assert len(result.rdatas) == 3
    assert res.shape == (12,)
    assert sres.shape == (12, problem.n_estimated)

    # least-squares identities (valid for parameter-independent sigmas)
    np.testing.assert_allclose(
        result.s2llh, sres.T @ sres, rtol=1e-5, atol=1e-8
    )
    np.testing.assert_allclose(
        [result.sllh[par_id] for par_id in x_ids],
        -res @ sres,
        rtol=1e-5,
        atol=1e-8,
    )

    # compare to finite differences
    ps.solver.set_sensitivity_order(SensitivityOrder.none)
    sres_fd = np.zeros_like(sres)
    for i, par_id in enumerate(x_ids):
        eps = 1e-5 * x[par_id]
        sres_fd[:, i] = (
            ps.simulate(x | {par_id: x[par_id] + eps}).res
            - ps.simulate(x | {par_id: x[par_id] - eps}).res
        ) / (2 * eps)
    np.testing.assert_allclose(sres, sres_fd, rtol=1e-4, atol=1e-6)

    # without sensitivities, no residual sensitivities are computed
    assert ps.simulate(x).sres is None


def test_aggregated_residual_sensitivities_reporting_modes(
    residual_test_importer,
):
    """Aggregated residuals and their sensitivities are available in the
    residual-only reporting mode, but not in the likelihood-only mode."""
    import numpy as np
    from amici.sim.sundials import (
        RDataReporting,
        SensitivityMethod,
        SensitivityOrder,
    )

    ps = residual_test_importer.create_simulator()
    ps.solver.set_sensitivity_method(SensitivityMethod.forward)
    ps.solver.set_sensitivity_order(SensitivityOrder.first)
    x = {"k_decay": 0.4, "scale": 1.5}

    expected = ps.simulate(x)

    ps.solver.set_return_data_reporting_mode(RDataReporting.residuals)
    result = ps.simulate(x)
    np.testing.assert_allclose(result.res, expected.res)
    np.testing.assert_allclose(result.sres, expected.sres)
    # the log-likelihood is computed from the residuals in this mode, ...
    np.testing.assert_allclose(result.llh, expected.llh)
    # ... but its sensitivities are not
    assert result.sllh is None
    assert result.s2llh is None

    ps.solver.set_return_data_reporting_mode(RDataReporting.likelihood)
    result = ps.simulate(x)
    assert result.res is None
    assert result.sres is None
    assert result.sllh is not None


def test_aggregated_residuals_failed_simulation(residual_test_importer):
    """No residuals are reported if any experiment failed to simulate.

    AMICI does not invalidate the residuals of failed simulations, so the
    residuals of the failed timepoints would look like a perfect fit.
    """
    import numpy as np
    from amici.sim.sundials import AMICI_SUCCESS, SensitivityMethod

    ps = residual_test_importer.create_simulator()
    ps.solver.set_sensitivity_method(SensitivityMethod.forward)
    ps.solver.set_sensitivity_order(SensitivityOrder.first)
    # provoke an integration failure
    ps.solver.set_max_steps(1)

    result = ps.simulate({"k_decay": 0.4, "scale": 1.5})
    assert any(rdata.status != AMICI_SUCCESS for rdata in result.rdatas)
    assert np.isnan(result.llh)
    assert result.res is None
    assert result.sres is None
    assert result.sllh is None
    assert result.s2llh is None


def test_aggregated_residual_sensitivities_non_gaussian_noise():
    """No residuals are computed for non-Gaussian noise models, but the
    log-likelihood and its sensitivities are."""
    import numpy as np
    from amici.sim.sundials import SensitivityMethod

    problem = Problem()
    problem.model = SbmlModel.from_antimony("xa = 1; xa' = -k_decay * xa;")
    problem.add_parameter(
        "k_decay", nominal_value=0.5, estimate=True, lb=0.01, ub=10
    )
    problem.add_observable(
        "obs_a",
        formula="xa",
        noise_formula="0.5",
        noise_distribution="laplace",
    )
    for t in (0.0, 1.0, 2.0):
        problem.add_measurement("obs_a", time=t, measurement=0.9 - 0.1 * t)
    problem.assert_valid()

    with TemporaryDirectoryWinSafe(prefix="petab_v2_laplace_") as output_dir:
        ps = PetabImporter(
            problem,
            module_name="test_petab_v2_residuals_laplace",
            output_dir=output_dir,
            verbose=False,
        ).create_simulator(force_import=True)
        ps.solver.set_sensitivity_method(SensitivityMethod.forward)
        ps.solver.set_sensitivity_order(SensitivityOrder.first)

        result = ps.simulate({"k_decay": 0.4})

    assert result.res is None
    assert result.sres is None
    # ... but the likelihood and its sensitivities are computed
    assert np.isfinite(result.llh)
    assert result.sllh.keys() == {"k_decay"}
    # the FIM is not available for non-Gaussian noise models
    assert result.s2llh is None


def test_jax_matches_sundials_with_per_observable_noise_parameters():
    """JAX and SUNDIALS agree on the log-likelihood for a model with
    per-observable noise placeholders (``noiseParameter1_{observableId}``).

    Regression test for a bug where the JAX backend applied the first
    observable's noise value to every measurement, since the PEtab-standard
    observable-suffixed placeholder names were not recognized as such and
    leaked into the model parameters instead of the per-measurement noise
    parameter array.
    """
    import numpy as np
    from petab.v2.core import ProblemConfig
    from petab.v2.models.sbml_model import SbmlModel

    problem = Problem()
    problem.config = ProblemConfig()
    problem.model = SbmlModel.from_antimony(
        "xa = 1; xb = 2; xa' = -0.1*xa; xb' = 0.2*xb;"
    )
    problem.add_parameter(
        "sd_obsA", nominal_value=0.5, estimate=True, lb=0.01, ub=10
    )
    problem.add_parameter(
        "sd_obsB", nominal_value=2.0, estimate=True, lb=0.01, ub=10
    )
    problem.add_observable(
        "obsA",
        "xa",
        noise_formula="noiseParameter1_obsA",
        noise_placeholders=["noiseParameter1_obsA"],
    )
    problem.add_observable(
        "obsB",
        "xb",
        noise_formula="noiseParameter1_obsB",
        noise_placeholders=["noiseParameter1_obsB"],
    )
    for t in [0.0, 1.0, 2.0]:
        problem.add_measurement(
            "obsA",
            time=t,
            measurement=1.0 + 0.05 * t,
            noise_parameters="sd_obsA",
        )
        problem.add_measurement(
            "obsB",
            time=t,
            measurement=2.0 + 0.1 * t,
            noise_parameters="sd_obsB",
        )

    from amici.sim.jax.petab import run_simulations

    # give the model a unique name -- the default derives from the (here
    # anonymous) SBML model id, which would collide with other tests using
    # an unnamed antimony model in this file, causing a stale cached module
    # to be reused
    jax_problem = PetabImporter(
        problem,
        jax=True,
        module_name="test_noise_params_jax",
        verbose=False,
    ).create_simulator(force_import=True)
    llh_jax, _ = run_simulations(jax_problem, max_steps=10_000)
    llh_sundials = (
        PetabImporter(
            problem,
            jax=False,
            module_name="test_noise_params_sundials",
            verbose=False,
        )
        .create_simulator(force_import=True)
        .simulate()
        .llh
    )
    np.testing.assert_allclose(float(llh_jax), float(llh_sundials), rtol=1e-5)


def test_jax_matches_sundials_with_arbitrarily_named_noise_placeholders():
    """JAX and SUNDIALS agree on the log-likelihood when noise placeholders
    do not follow the ``noiseParameter{n}_{observableId}`` naming
    convention.

    PEtab v2 declares observable/noise placeholders explicitly via the
    ``observablePlaceholders``/``noisePlaceholders`` columns, matched to
    measurement overrides purely by count and position (see
    ``petab.v2.lint.CheckOverridesMatchPlaceholders``) -- any name is
    allowed, e.g. the PEtab test suite (v2.0.0/sbml/0021) uses
    ``obs_a_noise_scale`` rather than ``noiseParameter1_obs_a``. Regression
    test for a bug where placeholders were only recognized when they
    happened to follow the suffixed naming convention.
    """
    import numpy as np
    from petab.v2.core import ProblemConfig
    from petab.v2.models.sbml_model import SbmlModel

    problem = Problem()
    problem.config = ProblemConfig()
    problem.model = SbmlModel.from_antimony(
        "xa = 1; xb = 2; xa' = -0.1*xa; xb' = 0.05*xb;"
    )
    problem.add_parameter(
        "scale_val", nominal_value=2.0, estimate=True, lb=0.01, ub=10
    )
    problem.add_parameter(
        "sd_b", nominal_value=1.0, estimate=True, lb=0.01, ub=10
    )
    problem.add_observable(
        "obs_a",
        formula="xa",
        noise_formula="obs_a * obs_a_noise_scale",
        noise_placeholders=["obs_a_noise_scale"],
    )
    problem.add_observable(
        "obs_b",
        formula="xb",
        noise_formula="sd_b_placeholder",
        noise_placeholders=["sd_b_placeholder"],
    )
    for t in [0.0, 1.0, 2.0]:
        problem.add_measurement(
            "obs_a",
            time=t,
            measurement=1.0 - 0.05 * t,
            noise_parameters="scale_val",
        )
        problem.add_measurement(
            "obs_b",
            time=t,
            measurement=2.0 + 0.05 * t,
            noise_parameters="sd_b",
        )

    from amici.sim.jax.petab import run_simulations

    jax_problem = PetabImporter(
        problem,
        jax=True,
        module_name="test_arbitrary_placeholder_jax",
        verbose=False,
    ).create_simulator(force_import=True)
    llh_jax, _ = run_simulations(jax_problem, max_steps=10_000)
    llh_sundials = (
        PetabImporter(
            problem,
            jax=False,
            module_name="test_arbitrary_placeholder_sundials",
            verbose=False,
        )
        .create_simulator(force_import=True)
        .simulate()
        .llh
    )
    np.testing.assert_allclose(float(llh_jax), float(llh_sundials), rtol=1e-5)


def _parameter_scales_test_problem() -> Problem:
    """A small two-experiment problem for parameter scale tests.

    The scaling factor ``s`` and the noise parameter ``sigma`` enter via
    placeholders that are overridden by problem parameters in ``e1``, but by
    numbers in ``e2``. ``k1`` is an estimated model parameter, ``k2`` a
    non-estimated one, which is turned into a fixed model parameter.
    """
    problem = Problem()
    problem.model = SbmlModel.from_antimony(
        """
        A = 1; B = 0
        k1 = 0.5; k2 = 0.1
        R1: A -> B; k1 * A
        R2: B -> A; k2 * B
        """
    )
    problem.add_observable(
        "obs_B",
        "observableParameter1_obs_B * B",
        noise_formula="noiseParameter1_obs_B",
        observable_placeholders=["observableParameter1_obs_B"],
        noise_placeholders=["noiseParameter1_obs_B"],
    )
    problem.add_parameter(
        "k1", estimate=True, nominal_value=0.5, lb=1e-3, ub=1e3
    )
    problem.add_parameter(
        "s", estimate=True, nominal_value=2.0, lb=1e-3, ub=1e3
    )
    problem.add_parameter("k2", estimate=False, nominal_value=0.1)
    problem.add_parameter("sigma", estimate=False, nominal_value=0.1)
    problem.add_condition("c0", A=1.0)
    problem.add_experiment("e1", 0, "c0")
    problem.add_experiment("e2", 0, "c0")
    for t, y in ((1.0, 0.8), (2.0, 1.0), (5.0, 1.4)):
        problem.add_measurement(
            "obs_B",
            experiment_id="e1",
            time=t,
            measurement=y,
            observable_parameters=["s"],
            noise_parameters=["sigma"],
        )
        problem.add_measurement(
            "obs_B",
            experiment_id="e2",
            time=t,
            measurement=y / 2,
            observable_parameters=[1.0],
            noise_parameters=[0.2],
        )
    problem.assert_valid()
    return problem


@pytest.fixture(scope="module")
def parameter_scales_test_importer() -> PetabImporter:
    """Importer for the model of :func:`_parameter_scales_test_problem`."""
    with TemporaryDirectoryWinSafe(prefix="petab_v2_pscale_") as output_dir:
        pi = PetabImporter(
            _parameter_scales_test_problem(),
            module_name="test_petab_v2_pscale",
            output_dir=output_dir,
            verbose=False,
        )
        pi.import_module(force_import=True)
        yield pi


#: Parameter scales for the tests based on
#: :func:`_parameter_scales_test_problem`
_PSCALE_TEST_SCALES = {
    "k1": "log10",
    "s": "log",
    "k2": "log10",
    "sigma": "log10",
}
#: Linear parameter values for the tests based on
#: :func:`_parameter_scales_test_problem`
_PSCALE_TEST_X = {"k1": 0.4, "s": 1.5, "k2": 0.2, "sigma": 0.3}


@pytest.fixture
def parameter_scales_simulators(parameter_scales_test_importer):
    """Simulators for :func:`_parameter_scales_test_problem` on linear scale
    and on the scales :data:`_PSCALE_TEST_SCALES`, with forward
    sensitivities and tight tolerances."""
    from amici.sim.sundials import SensitivityMethod

    simulators = (
        parameter_scales_test_importer.create_simulator(),
        parameter_scales_test_importer.create_simulator(
            parameter_scales=_PSCALE_TEST_SCALES
        ),
    )
    for ps in simulators:
        ps.solver.set_sensitivity_method(SensitivityMethod.forward)
        ps.solver.set_sensitivity_order(SensitivityOrder.first)
        ps.solver.set_absolute_tolerance(1e-14)
        ps.solver.set_relative_tolerance(1e-12)
    return simulators


def test_parameter_scales_default_to_linear(parameter_scales_test_importer):
    """Parameters without an explicitly chosen scale are on linear scale."""
    em = parameter_scales_test_importer.create_simulator().exp_man
    assert em.parameter_scales == dict.fromkeys(_PSCALE_TEST_SCALES, "lin")

    em = parameter_scales_test_importer.create_simulator(
        parameter_scales={"s": "log10"}
    ).exp_man
    assert em.parameter_scales == {
        "k1": "lin",
        "s": "log10",
        "k2": "lin",
        "sigma": "lin",
    }


def test_scale_unscale_parameters(parameter_scales_test_importer):
    """Parameter values are converted between linear scale and the chosen
    scales."""
    import numpy as np

    em = parameter_scales_test_importer.create_simulator(
        parameter_scales=_PSCALE_TEST_SCALES
    ).exp_man
    assert em.parameter_scales == _PSCALE_TEST_SCALES

    x_scaled = em.scale_parameters(_PSCALE_TEST_X)
    assert x_scaled == pytest.approx(
        {
            "k1": np.log10(0.4),
            "s": np.log(1.5),
            "k2": np.log10(0.2),
            "sigma": np.log10(0.3),
        }
    )
    assert em.unscale_parameters(x_scaled) == pytest.approx(_PSCALE_TEST_X)
    # integer-valued inputs, including negative exponents
    assert em.unscale_parameters(
        {"k1": -1, "s": 0, "k2": -2}
    ) == pytest.approx({"k1": 0.1, "s": 1.0, "k2": 0.01})


def test_parameter_scales_apply_parameters(parameter_scales_test_importer):
    """Each model parameter takes the scale of the problem parameter it is
    mapped to in the respective experiment. Everything else, including the
    fixed parameters, is on linear scale."""
    from amici.sim.sundials import ParameterScaling

    em = parameter_scales_test_importer.create_simulator(
        parameter_scales=_PSCALE_TEST_SCALES
    ).exp_man
    model = em.model
    pscales = {}
    for edata in em.create_edatas():
        em.apply_parameters(edata, em.scale_parameters(_PSCALE_TEST_X))
        pscales[edata.id] = dict(
            zip(model.get_free_parameter_ids(), edata.pscale)
        )
        fixed_parameters = dict(
            zip(model.get_fixed_parameter_ids(), edata.fixed_parameters)
        )
        assert fixed_parameters["k2"] == pytest.approx(_PSCALE_TEST_X["k2"])

    # `s` and `sigma` override the placeholders in `e1`, numbers in `e2`
    assert pscales == {
        "e1": {
            "k1": ParameterScaling.log10,
            "observableParameter1_obs_B": ParameterScaling.ln,
            "noiseParameter1_obs_B": ParameterScaling.log10,
        },
        "e2": {
            "k1": ParameterScaling.log10,
            "observableParameter1_obs_B": ParameterScaling.none,
            "noiseParameter1_obs_B": ParameterScaling.none,
        },
    }


def test_parameter_scales_sensitivities(parameter_scales_simulators):
    """Simulation results on non-linear parameter scales equal those on
    linear scale, and their sensitivities are related by the chain rule."""
    import numpy as np

    ps_lin, ps = parameter_scales_simulators
    x = _PSCALE_TEST_X
    result_lin = ps_lin.simulate(x)
    result = ps.simulate(ps.exp_man.scale_parameters(x))

    np.testing.assert_allclose(result.llh, result_lin.llh, rtol=1e-12)
    np.testing.assert_allclose(result.res, result_lin.res, rtol=1e-12)

    x_free_ids = ps.exp_man.petab_problem.x_free_ids
    assert x_free_ids == ["k1", "s"]
    # d x_lin / d x_scaled
    dlin_dscaled = np.array([x["k1"] * np.log(10), x["s"]])
    np.testing.assert_allclose(
        [result.sllh[par_id] for par_id in x_free_ids],
        [result_lin.sllh[par_id] for par_id in x_free_ids] * dlin_dscaled,
        rtol=1e-8,
    )
    np.testing.assert_allclose(
        result.sres, result_lin.sres * dlin_dscaled, rtol=1e-8, atol=1e-12
    )
    np.testing.assert_allclose(
        result.s2llh,
        result_lin.s2llh * np.outer(dlin_dscaled, dlin_dscaled),
        rtol=1e-8,
    )


def test_parameter_scales_simulate(parameter_scales_simulators):
    """Parameter values passed to the simulator are on the chosen scales,
    and unspecified parameters take their nominal values on those scales."""
    import numpy as np

    ps_lin, ps = parameter_scales_simulators
    np.testing.assert_allclose(
        ps.simulate().llh, ps_lin.simulate().llh, rtol=1e-12
    )
    np.testing.assert_allclose(
        ps.simulate({"k1": np.log10(0.4)}).llh,
        ps_lin.simulate({"k1": 0.4}).llh,
        rtol=1e-12,
    )
    # integer-valued inputs, including negative exponents
    np.testing.assert_allclose(
        ps.simulate({"k1": -1, "s": 0, "k2": -1, "sigma": -1}).llh,
        ps_lin.simulate({"k1": 0.1, "s": 1.0, "k2": 0.1, "sigma": 0.1}).llh,
        rtol=1e-12,
    )


def test_parameter_scales_invalid(parameter_scales_test_importer):
    """Invalid parameter scales are rejected."""
    with pytest.raises(ValueError, match="not PEtab problem parameters"):
        parameter_scales_test_importer.create_simulator(
            parameter_scales={
                "k1": "log10",
                "observableParameter1_obs_B": "log",
            }
        )
    with pytest.raises(ValueError, match="Invalid scale"):
        parameter_scales_test_importer.create_simulator(
            parameter_scales={"k1": "ln"}
        )


def test_parameter_scales_invalid_jax():
    """Invalid parameter scales are rejected for JAX models, too, before
    the model is imported."""
    with TemporaryDirectoryWinSafe(prefix="petab_v2_pscale_jax_") as tmp:
        pi = PetabImporter(
            _parameter_scales_test_problem(),
            module_name="test_petab_v2_pscale_jax",
            output_dir=tmp,
            jax=True,
            verbose=False,
        )
        with pytest.raises(ValueError, match="not PEtab problem parameters"):
            pi.create_simulator(parameter_scales={"typo": "lin"})
        with pytest.raises(ValueError, match="Invalid scale"):
            pi.create_simulator(parameter_scales={"k1": "ln"})
        with pytest.raises(NotImplementedError):
            pi.create_simulator(parameter_scales={"k1": "log10"})
