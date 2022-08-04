import amici
from util import create_sbml_model, create_amici_model
from test_pregenerated_models import (
    options_file, expected_results, expected_results_file,
    verify_simulation_results
)


def model_neuron_def():
    """Python implementation of the neuron model (Hodgkin-Huxley).

        ODEs
        ----
        d/dt v:
            - 0.04*v^2 + 5*v + 140 - u + I
        d/dt u:
            - a*(b*v - u);

        Events:
        -------
        event_1:
            trigger: v - 30
            bolus: [[ -c - v ],
                    [      0]]
        """
    # Model components
    species = ['v', 'u']
    initial_assignments = {
        'v': 'v0',
        'u': 'b*v0',
    }
    rate_rules = {
        'v': '0.04*v^2 + 5*v + 140 - u + I0',
        'u': 'a*(b*v - u)',
    }
    parameters = {
        'a': 0.02,
        'b': 0.3,
        'c': 65,
        'd': 0.9,
        'v0': -60,
        'I0': 10,
    }
    events = {
        'event_1': {
            'trigger': 'v > 30',
            'target': ['v', 'u'],
            'assignment': ['-c', 'd+u']
        },
    }

    observables = {
        'y1': {
            'name': 'v',
            'formula': 'v',
        }
    }

    event_observables = {
        'z1': {
            'name': 'z1',
            'event': 'event_1',
            'formula': 'time'
        }
    }
    return (
        initial_assignments,
        parameters,
        rate_rules,
        species,
        events,
        observables,
        event_observables
    )


def test_model_neuron():
    (
        initial_assignments,
        parameters,
        rate_rules,
        species,
        events,
        observables,
        event_observables
    ) = model_neuron_def()

    sbml_document, sbml_model = create_sbml_model(
        initial_assignments=initial_assignments,
        parameters=parameters,
        rate_rules=rate_rules,
        species=species,
        events=events,
        # uncomment `to_file` to save SBML model to file for inspection
        # to_file=sbml_test_models / (model_name + '.sbml'),
    )

    model = create_amici_model(
        sbml_model,
        model_name='model_neuron',
        observables=observables,
        constant_parameters=['v0', 'I0'],
        event_observables=event_observables
    )

    run_test_cases(model)

    return


def run_test_cases(model):

    solver = model.getSolver()

    model_name = model.getName()

    for case in list(expected_results[model_name].keys()):

        if case.startswith('sensi2'):
            continue

        amici.readModelDataFromHDF5(
            options_file, model.get(),
            f'/{model_name}/{case}/options'
        )
        amici.readSolverSettingsFromHDF5(
            options_file, solver.get(),
            f'/{model_name}/{case}/options'
        )

        edata = None
        if 'data' in expected_results[model.getName()][case].keys():
            edata = amici.readSimulationExpData(
                str(expected_results_file),
                f'/{model_name}/{case}/data', model.get()
            )
        rdata = amici.runAmiciSimulation(model, solver, edata)

        verify_simulation_opts = dict()

        if model_name.startswith('model_neuron'):
            verify_simulation_opts['atol'] = 1e-5
            verify_simulation_opts['rtol'] = 1e-2

        verify_simulation_results(
            rdata, expected_results[model.getName()][case]['results'],
            **verify_simulation_opts
        )


