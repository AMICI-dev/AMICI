"""Interface to facilitate AMICI generated models using JAX"""

from amici.jax.petab import JAXProblem, run_simulations
from amici.jax.model import JAXModel

__all__ = ["JAXModel", "JAXProblem", "run_simulations"]
