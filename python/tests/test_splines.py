import os
import numpy as np
import sympy as sp
import pandas as pd

from scipy.integrate import quad

import libsbml
import petab

from amici.sbml_utils import amici_time_symbol, setSbmlMath
from amici.splines import CubicHermiteSpline, UniformGrid


def create_spline_test_petab(yy_true, spline_dt=2.5, measure_dt=0.5, sigma=1.0, T=None, extrapolate=None, bc=None, folder=None):
    spline_dt = sp.nsimplify(spline_dt)
    tmax = (len(yy_true) - 1)*spline_dt
    xx = UniformGrid(0, tmax, spline_dt)
    yy = list(sp.symbols(f'y0:{len(yy_true)}'))

    spline = CubicHermiteSpline(
        'y', amici_time_symbol, xx, yy,
        bc=bc, extrapolate=extrapolate
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
    setSbmlMath(rule, spline.sbmlId)

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
    zdot = spline.formula.subs(dict(zip(yy, yy_true)))
    zdot = sp.lambdify(amici_time_symbol, zdot)
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

    if folder is not None:
        folder = os.path.abspath(folder)
        os.makedirs(folder, exist_ok=True)
        problem.to_files(
            sbml_file=os.path.join(folder, 'spline_test.xml'),
            condition_file=os.path.join(folder, 'spline_test_conditions.tsv'),
            measurement_file=os.path.join(folder, 'spline_test_measurements.tsv'),
            parameter_file=os.path.join(folder, 'spline_test_parameters.tsv'),
            observable_file=os.path.join(folder, 'spline_test_observables.tsv'),
            yaml_file=os.path.join(folder, 'spline_test.yaml')
        )

    return problem

if __name__ == "__main__":
    import sys
    folder = sys.argv[1] if len(sys.argv) > 1 else '.'
    yy_true = [0.0, 2.0, 3.0, 4.0, 1.0, -0.5, -1, -1.5, 0.5, 0.0]
    create_spline_test_petab(yy_true, folder=folder)
