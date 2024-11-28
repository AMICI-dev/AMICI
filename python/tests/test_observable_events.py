import os

import amici
import pytest
from test_pregenerated_models import (
    expected_results,
    expected_results_file,
    options_file,
    verify_simulation_results,
)
from util import create_amici_model, create_sbml_model


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
        observable: t
    """
    # Model components
    species = ["v", "u"]
    initial_assignments = {
        "v": "v0",
        "u": "b*v0",
    }
    rate_rules = {
        "v": "0.04*v^2 + 5*v + 140 - u + I0",
        "u": "a*(b*v - u)",
    }
    parameters = {
        "a": 0.02,
        "b": 0.3,
        "c": 65,
        "d": 0.9,
        "v0": -60,
        "I0": 10,
    }
    events = {
        "event_1": {
            "trigger": "v > 30",
            "target": ["v", "u"],
            "assignment": ["-c", "d+u"],
        },
    }

    observables = {
        "y1": {
            "name": "v",
            "formula": "v",
        }
    }

    event_observables = {
        "z1": {"name": "z1", "event": "event_1", "formula": "time"}
    }
    return (
        initial_assignments,
        parameters,
        rate_rules,
        species,
        events,
        observables,
        event_observables,
    )


def model_events_def():
    """Python implementation of the events model.

    ODEs
    ----
    d/dt x1:
        - -p1*heaviside(t-p4)*x1
    d/dt x2:
        - p2*x1*exp(-0.1*t)-p3*x2
    d/dt x3:
        - -x3+heaviside(t-4)

    Events:
    -------
    event_1:
        trigger: x2 > x3
        bolus: 0
        observable: t
    event_2:
        trigger: x1 > x3
        bolus: 0
        observable: t
    """
    # Model components
    species = ["x1", "x2", "x3"]
    initial_assignments = {
        "x1": "k1",
        "x2": "k2",
        "x3": "k3",
    }
    rate_rules = {
        "x1": "-p1*piecewise(1.0, time>p4, 0.0)*x1",
        "x2": "p2*x1*exp(-0.1*time)-p3*x2",
        "x3": "-x3+piecewise(1.0, time>4, 0.0)",
    }
    parameters = {
        "p1": 0.5,
        "p2": 2,
        "p3": 0.5,
        "p4": 0.5,
        "k1": 4,
        "k2": 8,
        "k3": 10,
        "k4": 4,
    }
    events = {
        "event_1": {"trigger": "x2 > x3", "target": [], "assignment": []},
        "event_2": {"trigger": "x1 > x3", "target": [], "assignment": []},
    }

    observables = {
        "y1": {
            "name": "y1",
            "formula": "p4*(x1+x2+x3)",
        }
    }

    event_observables = {
        "z1": {"name": "z1", "event": "event_1", "formula": "time"},
        "z2": {"name": "z2", "event": "event_2", "formula": "time"},
    }
    return (
        initial_assignments,
        parameters,
        rate_rules,
        species,
        events,
        observables,
        event_observables,
    )


models = [
    (model_neuron_def, "model_neuron_py", ["v0", "I0"]),
    (model_events_def, "model_events_py", ["k1", "k2", "k3", "k4"]),
]


@pytest.mark.skipif(
    os.environ.get("AMICI_SKIP_CMAKE_TESTS", "") == "TRUE",
    reason="skipping cmake based test",
)
@pytest.mark.parametrize("model_def,model_name,constants", models)
def test_models(model_def, model_name, constants):
    (
        initial_assignments,
        parameters,
        rate_rules,
        species,
        events,
        observables,
        event_observables,
    ) = model_def()

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
        model_name=model_name,
        observables=observables,
        constant_parameters=constants,
        event_observables=event_observables,
    )

    run_test_cases(model)

    return


def run_test_cases(model):
    solver = model.getSolver()

    model_name = model.getName()
    # we need a different name for the model module to avoid collisions
    #  with the matlab-pregenerated models, but we need the old name for
    #  the expected results
    model_name = model_name.removesuffix("_py")

    for case in list(expected_results[model_name].keys()):
        if case.startswith("sensi2"):
            continue

        amici.readModelDataFromHDF5(
            options_file, model.get(), f"/{model_name}/{case}/options"
        )
        amici.readSolverSettingsFromHDF5(
            options_file, solver.get(), f"/{model_name}/{case}/options"
        )

        edata = None
        if "data" in expected_results[model_name][case].keys():
            edata = amici.readSimulationExpData(
                str(expected_results_file),
                f"/{model_name}/{case}/data",
                model.get(),
            )
        rdata = amici.runAmiciSimulation(model, solver, edata)

        verify_simulation_opts = dict()

        if model_name.startswith("model_neuron"):
            verify_simulation_opts["atol"] = 1e-5
            verify_simulation_opts["rtol"] = 1e-2

        verify_simulation_results(
            rdata,
            expected_results[model_name][case]["results"],
            **verify_simulation_opts,
        )
