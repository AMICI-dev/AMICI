from amici.petab.petab_problem import AmiciPetabProblem
from benchmark_models_petab import get_problem


def test_amici_petab_problem_pregenerate():
    """AmiciPetabProblem with pre-generated ExpDatas"""
    petab_problem = get_problem("Boehm_JProteomeRes2014")

    # generate all edatas upon construction
    app = AmiciPetabProblem(petab_problem, store_edatas=True)
    assert len(app._edatas) == len(
        petab_problem.get_simulation_conditions_from_measurement_df()
    )

    for i, (_, condition) in enumerate(
        petab_problem.get_simulation_conditions_from_measurement_df().iterrows()
    ):
        assert app.get_edata(condition.simulationConditionId) is app._edatas[i]

    # ensure parameter are updated
    app.set_parameters(
        {app.model.getParameterIds()[0]: 0.12345}, scaled_parameters=True
    )
    for edata in app._edatas:
        assert edata.parameters[0] == 0.12345


def test_amici_petab_problem_on_demand():
    """AmiciPetabProblem with on-demand ExpDatas"""
    petab_problem = get_problem("Boehm_JProteomeRes2014")

    # generate all edatas upon construction
    app = AmiciPetabProblem(petab_problem, store_edatas=False)
    assert not app._edatas

    assert len(app.get_edatas()) == len(
        petab_problem.get_simulation_conditions_from_measurement_df()
    )

    # ensure parameter are updated
    app.set_parameters(
        {app.model.getParameterIds()[0]: 0.12345}, scaled_parameters=True
    )
    for edata in app.get_edatas():
        assert edata.parameters[0] == 0.12345

    some_sim_condition = (
        petab_problem.measurement_df.simulationConditionId.iloc[0]
    )
    # different objects
    assert app.get_edata(some_sim_condition) is not app.get_edata(
        some_sim_condition
    )
