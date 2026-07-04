import copy

from amici.importers.petab import (
    PetabImporter,
    flatten_timepoint_specific_output_overrides,
    has_timepoint_specific_overrides,
    unflatten_simulation_df,
)
from amici.sim.sundials import SensitivityOrder
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
    problem.model = SbmlModel.from_antimony("xx = 1; xx' = kk;")
    problem.add_parameter("kk", nominal_value=1.0, estimate=True, lb=0, ub=10)
    problem.add_observable("obs1", "xx", noise_formula="1")
    for i in range(5):
        problem.add_measurement("obs1", time=i, measurement=2 * i)

    pi = PetabImporter(problem)
    ps = pi.create_simulator(force_import=False)
    ps.solver.set_sensitivity_order(SensitivityOrder.none)

    ps_copy = copy.deepcopy(ps)

    assert ps.simulate({"kk": 2}).llh == ps_copy.simulate({"kk": 2}).llh

    ps.solver.set_sensitivity_order(SensitivityOrder.first)
    assert (
        ps.solver.get_sensitivity_order()
        != ps_copy.solver.get_sensitivity_order()
    )

    import pickle

    ps_pickle = pickle.loads(pickle.dumps(ps))
    assert ps.simulate({"kk": 2}).llh == ps_pickle.simulate({"kk": 2}).llh


def test_petab_importer_upgrades_v1_problem():
    """A PEtab v1 problem passed to PetabImporter is upgraded to PEtab v2."""
    import numpy as np
    import pandas as pd
    import petab.v1 as petab_v1
    import petab.v2
    from petab.v1.C import CONDITION_ID as V1_CONDITION_ID
    from petab.v1.C import (
        ESTIMATE,
        LIN,
        LOWER_BOUND,
        MEASUREMENT,
        NOISE_FORMULA,
        NOMINAL_VALUE,
        OBSERVABLE_FORMULA,
        OBSERVABLE_ID,
        PARAMETER_ID,
        PARAMETER_SCALE,
        SIMULATION_CONDITION_ID,
        TIME,
        UPPER_BOUND,
    )
    from petab.v1.models.sbml_model import SbmlModel as SbmlModelV1

    v1_problem = petab_v1.Problem(
        model=SbmlModelV1.from_antimony("xx = 1; xx' = kk;"),
        observable_df=pd.DataFrame(
            {
                OBSERVABLE_ID: ["obs1"],
                OBSERVABLE_FORMULA: ["xx"],
                NOISE_FORMULA: [1.0],
            }
        ).set_index(OBSERVABLE_ID),
        measurement_df=pd.DataFrame(
            {
                OBSERVABLE_ID: ["obs1"] * 3,
                SIMULATION_CONDITION_ID: ["c0"] * 3,
                TIME: [0.0, 1.0, 2.0],
                MEASUREMENT: [1.0, 2.0, 3.0],
            }
        ),
        condition_df=pd.DataFrame({V1_CONDITION_ID: ["c0"]}).set_index(
            V1_CONDITION_ID
        ),
        parameter_df=pd.DataFrame(
            {
                PARAMETER_ID: ["kk"],
                PARAMETER_SCALE: [LIN],
                LOWER_BOUND: [0.01],
                UPPER_BOUND: [10.0],
                NOMINAL_VALUE: [1.0],
                ESTIMATE: [1],
            }
        ).set_index(PARAMETER_ID),
    )

    importer = PetabImporter(v1_problem, verbose=False)
    # the in-memory v1 problem is upgraded to v2
    assert isinstance(importer.petab_problem, petab.v2.Problem)

    simulator = importer.create_simulator(force_import=True)
    assert np.isfinite(simulator.simulate().llh)


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
