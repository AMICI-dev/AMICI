"""Tests for petab_simulate.py."""

import os
import pytest
import tempfile

from amici.petab_simulate import PetabSimulator
import petab


@pytest.fixture
def petab_problem() -> petab.Problem:
    """Create a PEtab problem for use in tests."""
    petab_test_suite_path = os.path.join('..', '..', 'petab_test_suite')
    case = '0001'
    case_yaml_path = os.path.join('cases', case, f'_{case}.yaml')
    petab_yaml_path = os.path.join(petab_test_suite_path,
                                   case_yaml_path)
    return petab.Problem.from_yaml(petab_yaml_path)


def test_simulate_without_noise(petab_problem):
    """Test the reproducibility of simulation without noise."""
    simulator = PetabSimulator(petab_problem)
    print(simulator.working_dir)
    synthetic_data_df_a = simulator.simulate()
    synthetic_data_df_b = simulator.simulate()
    simulator.remove_working_dir()
    # Repeated simulation without noise is reproducible.
    assert synthetic_data_df_b.equals(synthetic_data_df_a)

    simulator = PetabSimulator(petab_problem)
    print(simulator.working_dir)
    synthetic_data_df_c = simulator.simulate()
    simulator.remove_working_dir()
    # Repeated simulation without noise, with a different PetabSimulator
    # instance, is reproducible.
    assert synthetic_data_df_c.equals(synthetic_data_df_a)


def test_subset_call(petab_problem):
    """
    Test the ability to customize AMICI methods, specifically:
    :py:func:`amici.petab_import.import_petab_problem` (`model_name`,
    `model_output_dir`, import is skipped if `amici_model` is specified), and
    :py:func:`amici.petab_objective.simulate_petab` (`amici_model`, `solver`).
    """
    model_name = 'model_name_dummy'
    model_output_dir = tempfile.mkdtemp()

    simulator0 = PetabSimulator(petab_problem)
    assert not os.path.isdir(os.path.join(model_output_dir, model_name))
    simulator0.simulate(model_name=model_name,
                        model_output_dir=model_output_dir)
    # Model name is handled correctly
    assert simulator0.amici_model.getName() == model_name
    # Check model output directory is created, by
    # :py:func:`amici.petab_import.import_petab_problem`
    assert os.path.isdir(os.path.join(model_output_dir, model_name))

    simulator = PetabSimulator(petab_problem)
    simulator.simulate(amici_model=simulator0.amici_model)
    # AMICI model is handled correctly to skip import
    assert simulator.amici_model.getName() == model_name
    simulator.simulate()
    # AMICI model persists between :py:func:`PetabSimulator.simulate` calls
    assert simulator.amici_model.getName() == model_name
    # Erroneous solver raises an error
    with pytest.raises(TypeError):
        simulator.simulate(solver=False)

    simulator0.remove_working_dir()
    simulator.remove_working_dir()
