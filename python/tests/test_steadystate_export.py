"""Tests for steadystate export"""

import petab
import petabtests
import amici
import amici.petab_import
from amici.petab_objective import export_steadystates
import numpy as np
import pandas as pd
import warnings
import os


def create_petab_model():
    """Get a PEtab problem for testing the steady state export"""

    # Use case 0009 from the PEtab test suite as test model
    case = petabtests.test_id_str('0009')
    case_dir = os.path.join(os.path.dirname(petabtests.__file__), '..',
                            'cases', case)
    # import PEtab problem
    yaml_file = os.path.join(case_dir, petabtests.problem_yaml_name(case))
    petab_problem = petab.Problem.from_yaml(yaml_file)

    # compile amici model
    model_output_dir = 'amici_models/model_0009'
    amici_model = amici.petab_import.import_petab_problem(
        petab_problem, model_output_dir=model_output_dir,
        force_compile=True)

    return petab_problem, amici_model

def extend_petab_model(petab_problem):
    """Add rows with steady state measurements to the tables"""
    # apppend a condition
    c_add = pd.DataFrame(columns=('k1',), index=('c1',), data=np.array((1.,)))
    petab_problem.condition_df = petab_problem.condition_df.append(c_add)

    # append four measurements for two conditions
    m_data = petab_problem.measurement_df.loc[:].append(
        petab_problem.measurement_df.loc[:])
    m_data.index = index=range(2,6)
    m_data['preequilibrationConditionId'] = ['preeq_c0'] * 3 + ['']
    m_data['simulationConditionId'] = ['c1', 'c1', 'c1', 'preeq_c0']
    m_data['time'] = [5.0, float('inf'), float('inf'), float('inf')]
    m_data['measurement'] = [0.4, 0.9, 0.8, 0.2]
    m_add = pd.DataFrame(columns=petab_problem.measurement_df.columns,
                         data=m_data)
    petab_problem.measurement_df = petab_problem.measurement_df.append(m_add)

    # return the updated problem
    return petab_problem


def test_steadystate_export():
    """Exporting steady states"""

    # Get the PEtab problem first
    petab_problem, amici_model  = create_petab_model()
    
    # test passing nothing but a petab problem, having no steady state data
    steady_states = export_steadystates(petab_problem, amici_model)
    assert (steady_states.values.size == 0)

    # add steady state data
    petab_problem = extend_petab_model(petab_problem)
    steady_states = export_steadystates(petab_problem, amici_model)
    assert (steady_states.values.shape == (2,2))
    assert [True, False] == [':' in cond for cond in steady_states.index]

    # test exporting the csv
    tmp_file = os.path.join(os.path.abspath('amici_models/model_0009/'),
                            'deleteme_ss.csv')
    assert not os.path.exists(tmp_file)
    _ = export_steadystates(petab_problem, amici_model, export_file=tmp_file)
    assert os.path.exists(tmp_file)
    os.remove(tmp_file)
    
    # test passing simulation conditions
    sim_cond = petab_problem.get_simulation_conditions_from_measurement_df()
    sim_cond_1 = sim_cond.drop(1)
    steady_states = export_steadystates(petab_problem, amici_model,
                                        simulation_conditions=sim_cond_1)
    assert (steady_states.values.shape == (1,2))
    assert ([False,] == [':' in cond for cond in steady_states.index])

    # test passing edatas
    sim_cond_2 = sim_cond.drop(2)
    edatas = amici.petab_objective.create_edatas(amici_model, petab_problem,
                                                 sim_cond_2)
    steady_states = export_steadystates(petab_problem, amici_model,
                                        simulation_conditions=sim_cond_2,
                                        edatas=edatas)
    assert (steady_states.values.shape == (1,2))
    assert (['edata_1',] == list(steady_states.index))

    # test passing only conditions without steadystate measurements
    assert not os.path.exists(tmp_file)
    sim_cond_0 = sim_cond.drop([1,2], axis=0)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        steady_states = export_steadystates(petab_problem, amici_model,
                                            simulation_conditions=sim_cond_0,
                                            export_file=tmp_file)
        assert (steady_states.values.size == 0)
        assert not os.path.exists(tmp_file)
