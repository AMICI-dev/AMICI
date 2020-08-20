"""Test support for special noise distributions."""
import numpy as np
from scipy.special import loggamma
import pandas as pd

from amici.petab_import import import_petab_problem
from amici.petab_objective import (
        simulate_petab, rdatas_to_measurement_df)
import petab


def test_binomial():
    """Test a model with binomial noise model."""
    problem = petab.Problem.from_yaml('binomial_model/_binomial_model.yaml')
    model = import_petab_problem(problem, force_compile=True)

    # measurements
    measurement_df = problem.measurement_df

    # simulate
    ret = simulate_petab(problem, model)
    
    rdatas = ret['rdatas']
    llh = ret['llh']
    simulation_df = rdatas_to_measurement_df(rdatas, model, measurement_df)

    llh_exp = binomial_llh(
        simulation_df.measurement, measurement_df.measurement, 0.5)

    assert np.isclose(llh, llh_exp)


def binomial_llh(simu, meas, p):
    """Compute a binomial llh."""
    if any(meas > simu):
        return - np.inf
    llhs = loggamma(simu+1) - loggamma(meas+1) - loggamma(simu-meas+1) \
        + meas * np.log(p) + (simu-meas) * np.log(1-p)
    return sum(llhs)
