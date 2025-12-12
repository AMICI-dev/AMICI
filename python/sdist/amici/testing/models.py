"""Test models used by AMICI."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

import libsbml
import sympy as sp

from amici import MeasurementChannel, SbmlImporter, import_model_module
from amici.importers.antimony import antimony2amici, antimony2sbml
from amici.importers.sbml.splines import CubicHermiteSpline
from amici.importers.utils import amici_time_symbol

if TYPE_CHECKING:
    from amici.sim.sundials import AmiciModel, Model

model_dirac_ant = r"""
p1 = 1;
p2 = 0.5;
p3 = 2;
p4 = 3;

x1 = 0;
x2 = 0;

x1' = -p1 * x1
x2' = p3 * x1 - p4 * x2;

at time >= p2: x1 = x1 + 1
"""


calvetti_ant = r"""
# constants
V1ss = 0.29;
R1ss = 0.74;
V2ss = 0.44;
R2ss = 0.08;
V3ss = 0.27;
R3ss = 0.18;

# static expressions
p1 = 1;
p2 = 1 - R1ss;
p3 = 1 - (R1ss + R2ss);
L1 = (R1ss * (V1ss ^ 2)) ^ (1 / 3);
L2 = (R2ss * (V2ss ^ 2)) ^ (1 / 3);
L3 = (R3ss * (V3ss ^ 2)) ^ (1 / 3);
C1ss = V1ss / (p1 - 0.5 * R1ss);
C2ss = V2ss / (p2 - 0.5 * R2ss);
C3ss = V3ss / (p3 - 0.5 * R3ss);

species V1, V2, V3, f1, f2, f3;
# initial values
V1 = V1ss;
V2 = V2ss;
V3 = V3ss;
f1 = 1;
f2 = 1;
f3 = 1;

R1 := (L1 ^ 3) / (V1 ^ 2);
R2 := (L2 ^ 3) / (V2 ^ 2);
R3 := (L3 ^ 3) / (V3 ^ 2);
f0 := (2 / R1) - ((R1 + R2) * f1 + (R2 + R3) * f2 + R3 * f3) / R1;
s := piecewise(0, time < 10, 1, time >= 10 & & time <= 12, 0);

# differential equations
rate_of_V1 := (1 / 31) * (((1.29 - (V1 / V1ss)) / (1.29 - 1)) + s - 2 * (V1 / (C1ss * ((R1 + R2) * f1 + (R2 + R3) * f2 + R3 * f3))));
V1' = rate_of_V1;
rate_of_V2 := (1 / 163) * (((1.51 - (V2 / V2ss)) / (1.51 - 1)) - 2 * (V2 / (C2ss * ((R2 + R3) * f2 + R3 * f3))));
V2' = rate_of_V2;
rate_of_V3 := (1 / 122) * (((1000000 - (V3 / V3ss)) / (1000000 - 1)) - 2 * (V3 / (C3ss * (R3 * f3))));
V3' = rate_of_V3;

# algebraic equations
# TODO: AMICI does currently not support rateOf(.) inside AlgebraicRules
# 0 = f0-rateOf(V1)-f1;
# 0 = f1-rateOf(V2)-f2;
# 0 = f2-rateOf(V3)-f3;
0 = f0 - rate_of_V1 - f1;
0 = f1 - rate_of_V2 - f2;
0 = f2 - rate_of_V3 - f3;
"""

robertson_ant = r"""
model robertson

p1 = 0.04;
p2 = 1e4;
p3 = 3e7;
k1 = 0.9;

x1 = k1;
x2 = 0;
x3 = 0;

# differential states
x1' = -p1 * x1 + p2 * x2 * x3;
x2' = p1 * x1 - p2 * x2 * x3 - p3 * x2^2;

# algebraic equations
0 = x1 + x2 + x3 - 1;

end
"""


def import_model_robertson(output_dir: Path = None) -> Model:
    """Import the Robertson model."""
    model_name = "model_robertson_py"

    if output_dir is None:
        output_dir = model_name

    antimony2amici(
        robertson_ant,
        fixed_parameters=["k1"],
        observation_model=[
            MeasurementChannel(id_="obs_x1", formula="x1"),
            MeasurementChannel(id_="obs_x2", formula="1e4 * x2"),
            MeasurementChannel(id_="obs_x3", formula="x3"),
        ],
        model_name=model_name,
        output_dir=output_dir,
    )
    model_module = import_model_module(model_name, output_dir)
    model = model_module.get_model()

    return model


def import_model_calvetti(output_dir: Path = None) -> Model:
    model_name = "model_calvetti_py"

    if output_dir is None:
        output_dir = model_name

    antimony2amici(
        calvetti_ant,
        fixed_parameters=["V1ss", "R1ss", "V2ss", "R2ss", "V3ss", "R3ss"],
        observation_model=[
            MeasurementChannel(id_="obs_V1", formula="V1"),
            MeasurementChannel(id_="obs_V2", formula="V2"),
            MeasurementChannel(id_="obs_V3", formula="V3"),
            MeasurementChannel(id_="obs_f0", formula="f0"),
            MeasurementChannel(id_="obs_f1", formula="f1"),
            MeasurementChannel(id_="obs_f2", formula="f2"),
        ],
        model_name=model_name,
        output_dir=output_dir,
        hardcode_symbols=["p1"],
    )
    model_module = import_model_module(model_name, output_dir)
    model = model_module.get_model()

    assert model.get_fixed_parameter_ids() == (
        "V1ss",
        "R1ss",
        "V2ss",
        "R2ss",
        "V3ss",
        "R3ss",
    )

    return model


def import_model_dirac(output_dir: Path = None) -> Model:
    """Import the Dirac model."""
    model_name = "model_dirac_py"

    if output_dir is None:
        output_dir = model_name

    antimony2amici(
        model_dirac_ant,
        observation_model=[MeasurementChannel(id_="obs_x2", formula="x2")],
        model_name=model_name,
        output_dir=output_dir,
    )
    model_module = import_model_module(model_name, output_dir)
    model = model_module.get_model()

    return model


def import_model_neuron(output_dir: Path = None) -> AmiciModel:
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

    sbml_document, sbml_model = create_sbml_model(
        initial_assignments=initial_assignments,
        parameters=parameters,
        rate_rules=rate_rules,
        species=species,
        events=events,
        # uncomment `to_file` to save SBML model to file for inspection
        # to_file=sbml_test_models / (model_name + '.sbml'),
    )

    model_name = "model_neuron_py"
    constants = ["v0", "I0"]
    model = create_amici_model(
        sbml_model,
        model_name=model_name,
        observation_model=[
            MeasurementChannel(id_="y1", name="v", formula="v"),
            MeasurementChannel(id_="z1", event_id="event_1", formula="time"),
        ],
        fixed_parameters=constants,
        output_dir=output_dir,
    )
    return model


def import_model_events(output_dir: Path = None) -> AmiciModel:
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

    sbml_document, sbml_model = create_sbml_model(
        initial_assignments=initial_assignments,
        parameters=parameters,
        rate_rules=rate_rules,
        species=species,
        events=events,
        # uncomment `to_file` to save SBML model to file for inspection
        # to_file=sbml_test_models / (model_name + '.sbml'),
    )

    model_name = "model_events_py"
    constants = ["k1", "k2", "k3", "k4"]

    model = create_amici_model(
        sbml_model,
        model_name=model_name,
        observation_model=[
            MeasurementChannel(id_="y1", formula="p4*(x1+x2+x3)", name="y1"),
            MeasurementChannel(id_="z1", event_id="event_1", formula="time"),
            MeasurementChannel(id_="z2", event_id="event_2", formula="time"),
        ],
        fixed_parameters=constants,
        output_dir=output_dir,
    )
    return model


def import_model_jakstat(output_dir: Path = None) -> AmiciModel:
    model_name = "model_jakstat_adjoint_py"
    if output_dir is None:
        output_dir = model_name

    # Create basic SBML model without spline
    ant_str = r"""
    model model_jakstat_adjoint
        # parameters
        p1 = 10^0.60;
        p2 = 10^3;
        p3 = 10^-0.95;
        p4 = 10^-0.0075;
        init_STAT = 10^0;
        # spline values
        sp1 = 10^-2.8;
        sp2 = 10^-0.26;
        sp3 = 10^-0.075;
        sp4 = 10^-0.41;
        sp5 = 10^-5;
        # output parameters
        offset_tSTAT = 10^-0.74;
        offset_pSTAT = 10^-0.64;
        scale_tSTAT = 10^-0.11;
        scale_pSTAT = 10^0.027;
        sigma_pSTAT = 10^-0.5;
        sigma_tSTAT = 10^0;
        sigma_pEpoR = 10^-0.5;

        # constants
        Omega_cyt = 1.4;
        Omega_nuc = 0.45;

        # state variables
        species STAT, pSTAT, pSTAT_pSTAT, npSTAT_npSTAT, nSTAT1, nSTAT2, nSTAT3, nSTAT4, nSTAT5;

        STAT = init_STAT;
        pSTAT = 0;
        pSTAT_pSTAT = 0;
        npSTAT_npSTAT = 0;
        nSTAT1 = 0;
        nSTAT2 = 0;
        nSTAT3 = 0;
        nSTAT4 = 0;
        nSTAT5 = 0;

        STAT' = (Omega_nuc*p4*nSTAT5 - Omega_cyt*STAT*p1*u)/Omega_cyt;
        pSTAT' = STAT*p1*u - 2*p2*pSTAT^2;
        pSTAT_pSTAT' = p2*pSTAT^2 - p3*pSTAT_pSTAT;
        npSTAT_npSTAT' = -(Omega_nuc*p4*npSTAT_npSTAT - Omega_cyt*p3*pSTAT_pSTAT)/Omega_nuc;
        nSTAT1' = -p4*(nSTAT1 - 2*npSTAT_npSTAT);
        nSTAT2' = p4*(nSTAT1 - nSTAT2);
        nSTAT3' = p4*(nSTAT2 - nSTAT3);
        nSTAT4' = p4*(nSTAT3 - nSTAT4);
        nSTAT5' = p4*(nSTAT4 - nSTAT5);

        # target of spline to be added below
        var u = 0;
    end
    """
    sbml_str = antimony2sbml(ant_str)
    sbml_doc = libsbml.SBMLReader().readSBMLFromString(sbml_str)
    sbml_model = sbml_doc.getModel()

    # Add spline to the SBML model
    # Timepoints for spline nodes
    nodes = [0, 5, 10, 20, 60]
    # Parameterized spline values at nodes
    values_at_nodes = [sp.Symbol(f"sp{i + 1}") for i, t in enumerate(nodes)]
    # The original MATLAB test used a natural cubic spline, but this is not
    # supported during Python import, so we use a cubic Hermite spline instead.
    # This is reasonably close to the original spline, but not exactly the
    # same.
    spline = CubicHermiteSpline(
        sbml_id="u",
        evaluate_at=amici_time_symbol,
        nodes=nodes,
        values_at_nodes=values_at_nodes,
        extrapolate=("constant", "constant"),
        logarithmic_parametrization=True,
    )
    spline.add_to_sbml_model(sbml_model, auto_add=False)

    SbmlImporter(sbml_model).sbml2amici(
        fixed_parameters=["Omega_cyt", "Omega_nuc"],
        observation_model=[
            MeasurementChannel(
                id_="obs_pSTAT",
                formula="offset_pSTAT + scale_pSTAT/init_STAT*(pSTAT + 2*pSTAT_pSTAT)",
                sigma="sigma_pSTAT",
            ),
            MeasurementChannel(
                id_="obs_tSTAT",
                formula="offset_tSTAT + scale_tSTAT/init_STAT*(STAT + pSTAT + 2*(pSTAT_pSTAT))",
                sigma="sigma_tSTAT",
            ),
            MeasurementChannel(
                id_="obs_spline", formula="u", sigma="sigma_pEpoR"
            ),
        ],
        model_name=model_name,
        output_dir=output_dir,
    )
    model_module = import_model_module(model_name, output_dir)
    model = model_module.get_model()
    return model


def import_model_nested_events(output_dir: Path = None) -> AmiciModel:
    model_name = "model_nested_events_py"
    if output_dir is None:
        output_dir = model_name

    ant_str = r"""
    model nested_events

    # parameters
    V_0 = 0.1;
    V_0_inject = 1000;
    t_0 = 2;
    rho_V = 8e-1;
    delta_V = 1.6;

    # state variables
    Virus = V_0
    Virus' = piecewise(0, Virus < 1, Virus*rho_V) -Virus*delta_V;

    injection: at time >= t_0: Virus = Virus + V_0_inject;

    end
    """

    antimony2amici(
        ant_str,
        observation_model=[
            MeasurementChannel(id_="obs_Virus", formula="Virus"),
        ],
        model_name=model_name,
        output_dir=output_dir,
    )
    model_module = import_model_module(model_name, output_dir)
    model = model_module.get_model()
    return model


def import_model_steadystate(output_dir: Path = None) -> AmiciModel:
    model_name = "model_steadystate_py"
    if output_dir is None:
        output_dir = model_name

    ant_str = r"""
    model model_steadystate

    k1 = 0.1;
    k2 = 0.4;
    k3 = 0.7;
    k4 = 1;

    p1 = 1;
    p2 = 0.5;
    p3 = 0.4;
    p4 = 2;
    p5 = 0.1;

    species x1 = k1;
    species x2 = k2;
    species x3 = k3;
    x1' = -2*p1*x1^2 - p2*x1*x2 + 2*p3*x2 + p4*x3 + p5;
    x2' = p1*x1^2 - p2*x1*x2 - p3*x2 + p4*x3;
    x3' = p2*x1*x2 - p4*x3 - k4*x3;

    end
    """

    antimony2amici(
        ant_str,
        fixed_parameters=["k1", "k2", "k3", "k4"],
        observation_model=[
            MeasurementChannel(id_="obs_x1", formula="x1"),
            MeasurementChannel(id_="obs_x2", formula="x2"),
            MeasurementChannel(id_="obs_x3", formula="x3"),
        ],
        model_name=model_name,
        output_dir=output_dir,
    )
    model_module = import_model_module(model_name, output_dir)
    model = model_module.get_model()
    return model


def import_test_models():
    """Import models required for C++ integration tests."""
    # out_root = Path(os.getcwd())
    repo_root = Path(__file__).parents[4]
    out_root = repo_root / "models"

    print(f"Generating test models in {out_root}...")
    print("Importing model_dirac_py...")
    import_model_dirac(output_dir=out_root / "model_dirac_py")
    print("Importing model_events_py...")
    import_model_events(output_dir=out_root / "model_events_py")
    print("Importing model_neuron_py...")
    import_model_neuron(output_dir=out_root / "model_neuron_py")
    print("Importing model_calvetti_py...")
    import_model_calvetti(output_dir=out_root / "model_calvetti_py")
    print("Importing model_robertson_py...")
    import_model_robertson(output_dir=out_root / "model_robertson_py")
    print("Importing model_jakstat_adjoint_py...")
    import_model_jakstat(output_dir=out_root / "model_jakstat_adjoint_py")
    print("Importing model_nested_events_py...")
    import_model_nested_events(output_dir=out_root / "model_nested_events_py")
    print("Importing model_steadystate_py...")
    import_model_steadystate(output_dir=out_root / "model_steadystate_py")


def create_sbml_model(
    initial_assignments,
    parameters,
    rate_rules,
    species,
    events,
    to_file: str = None,
):
    """Create an SBML model from simple definitions.

    See the model definitions and usage in :py:func:`model` for example input.

    The default initial concentration of species is `1.0`. This can currently
    be changed by specifying an initial assignment.
    """
    document = libsbml.SBMLDocument(3, 1)
    model = document.createModel()

    compartment = model.createCompartment()
    compartment.setId("compartment")
    compartment.setConstant(True)
    compartment.setSize(1)
    compartment.setSpatialDimensions(3)
    compartment.setUnits("dimensionless")

    for species_id in species:
        species = model.createSpecies()
        species.setId(species_id)
        species.setCompartment("compartment")
        species.setConstant(False)
        species.setSubstanceUnits("dimensionless")
        species.setBoundaryCondition(False)
        species.setHasOnlySubstanceUnits(False)
        species.setInitialConcentration(1.0)

    for target, formula in initial_assignments.items():
        initial_assignment = model.createInitialAssignment()
        initial_assignment.setSymbol(target)
        initial_assignment.setMath(libsbml.parseL3Formula(formula))

    for target, formula in rate_rules.items():
        rate_rule = model.createRateRule()
        rate_rule.setVariable(target)
        rate_rule.setMath(libsbml.parseL3Formula(formula))

    for parameter_id, parameter_value in parameters.items():
        parameter = model.createParameter()
        parameter.setId(parameter_id)
        parameter.setConstant(True)
        parameter.setValue(parameter_value)
        parameter.setUnits("dimensionless")

    for event_id, event_def in events.items():
        event = model.createEvent()
        event.setId(event_id)
        event.setName(event_id)
        event.setUseValuesFromTriggerTime(False)
        trigger = event.createTrigger()
        trigger.setMath(libsbml.parseL3Formula(event_def["trigger"]))
        trigger.setPersistent(True)
        trigger.setInitialValue(True)

        def create_event_assignment(target, assignment):
            ea = event.createEventAssignment()
            ea.setVariable(target)
            ea.setMath(libsbml.parseL3Formula(assignment))

        if isinstance(event_def["target"], list):
            for event_target, event_assignment in zip(
                event_def["target"], event_def["assignment"], strict=True
            ):
                create_event_assignment(event_target, event_assignment)

        else:
            create_event_assignment(
                event_def["target"], event_def["assignment"]
            )

    if to_file:
        libsbml.writeSBMLToFile(document, to_file)

    # Need to return document, else AMICI throws an error.
    # (possibly due to garbage collection?)
    return document, model


def create_amici_model(
    sbml_model, model_name, output_dir: Path = None, **kwargs
) -> AmiciModel:
    """
    Import an sbml file and create an AMICI model from it
    """
    sbml_importer = SbmlImporter(sbml_model)

    if output_dir is None:
        sbml_test_models_output_dir = Path("amici_models")
        sbml_test_models_output_dir.mkdir(parents=True, exist_ok=True)
        # try not to exceed the stupid maximum path length on windows 💩
        output_dir = (
            sbml_test_models_output_dir / model_name
            if sys.platform != "win32"
            else tempfile.mkdtemp()
        )

    sbml_importer.sbml2amici(
        model_name=model_name, output_dir=output_dir, **kwargs
    )

    model_module = import_model_module(model_name, output_dir)
    return model_module.get_model()
