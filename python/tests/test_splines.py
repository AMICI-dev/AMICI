import pandas as pd
from scipy.integrate import quad
import petab

def create_spline_test_petab(yy_true: Sequence, spline_dt=2.5, measure_dt=0.5, sigma=1.0, T=None, extrapolate=True, periodic=True):
    spline_dt = sp.nsimplify(spline_dt)
    tmax = (len(yy_true) - 1)*spline_dt
    yy = list(sp.symbols(f'y0:{len(yy_true)}'))
    spline = CatmullRomSpline_Piecewise(
        'y',
        AMICI_TIME,
        UniformGrid(0, tmax, spline_dt),
        yy,
        extrapolate=extrapolate,
        periodic=periodic
    )

    # Create SBML document
    doc = libsbml.SBMLDocument(2, 5)
    model = doc.createModel()
    model.setId('spline_test')

    # Add compartment
    compartmentId = 'compartment'
    comp = model.createCompartment()
    comp.setId(compartmentId)
    comp.setSize(1)

    # Add spline to the model
    spline.addToSbmlModel(model, auto_add=True, y_nominal=yy_true)

    # Add a species
    speciesId = 'z'
    species = model.createSpecies()
    species.setId(speciesId)
    species.setCompartment(compartmentId)
    species.setInitialAmount(0.0)

    # The derivative of the species is equal to the spline
    rule = model.createRateRule()
    rule.setIdAttribute(f'derivative_of_{speciesId}')
    rule.setVariable(speciesId)
    setMath(rule, spline.sbmlId)

    # Create condition table
    condition_df = pd.DataFrame({'conditionId' : ['condition1']})
    condition_df.set_index(['conditionId'], inplace=True)

    # Create parameter table
    parameter_df = pd.DataFrame({
        'parameterId' : [y.name for y in yy],
        'parameterScale' : 'lin',
        'lowerBound' : -10,
        'upperBound' :  10,
        'nominalValue' : 0,
        'estimate' : 1
    })
    parameter_df.set_index(['parameterId'], inplace=True)

    # Create observable table
    observable_df = pd.DataFrame({
        'observableId' : ['z_obs'],
        'observableFormula' : ['z'],
        'observableTransformation' : ['lin'],
        'noiseFormula' : [sigma],
        'noiseDistribution' : ['normal']
    })
    observable_df.set_index(['observableId'], inplace=True)

    # Evaluate spline
    if T is None:
        T = tmax
    T = sp.nsimplify(T)
    measure_dt= sp.nsimplify(measure_dt)
    n_obs = int(T / measure_dt) + 1
    T = (n_obs - 1) * measure_dt
    t_obs = np.asarray(UniformGrid(0, T, measure_dt), dtype=float)
    zdot = spline.formula(sbml=False).subs(dict(zip(yy, yy_true)))
    zdot = sp.lambdify(AMICI_TIME, zdot)
    z_true = np.concatenate((
        [0],
        np.cumsum([
            quad(zdot, float(i * measure_dt), float((i+1) * measure_dt))[0]
            for i in range(n_obs - 1)
        ])
    ))
    if len(z_true) != len(t_obs):
        print(len(z_true))
        print(len(t_obs))
    z_obs = z_true + np.random.randn(len(z_true))

    # Create measurement table
    measurement_df = pd.DataFrame({
        'observableId' : 'z_obs',
        'simulationConditionId' : 'condition1',
        'time' : t_obs,
        'measurement' : z_obs
    })

    problem = petab.Problem(
        sbml_document = doc,
        sbml_model = model,
        condition_df = condition_df,
        measurement_df = measurement_df,
        parameter_df = parameter_df,
        observable_df = observable_df
    )

    if petab.lint_problem(problem):
        raise Exception('PEtab lint failed')

    return problem

problem = create_spline_test_petab(
    [0.0, 2.0, 3.0, 4.0, 1.0, -0.5, -1, -1.5, 0.5, 0.0],
    spline_dt=2.5,
    measure_dt=0.5,
    sigma=1.0,
    extrapolate=True,
    periodic=False
)

import os
root = ''
os.makedirs(root, exist_ok=True)
problem.to_files(
    sbml_file=os.path.join(root, 'spline_test.xml'),
    condition_file=os.path.join(root, 'spline_test_conditions.tsv'),
    measurement_file=os.path.join(root, 'spline_test_measurements.tsv'),
    parameter_file=os.path.join(root, 'spline_test_parameters.tsv'),
    observable_file=os.path.join(root, 'spline_test_observables.tsv'),
    yaml_file=os.path.join(root, 'spline_test.yaml')
)
