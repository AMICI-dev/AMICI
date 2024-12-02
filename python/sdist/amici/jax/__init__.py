"""Interface to facilitate AMICI generated models using JAX"""

from amici.jax.petab import JAXProblem, run_simulations
from amici.jax.model import JAXModel
from amici.jax.nn import generate_equinox

__all__ = ["JAXModel", "JAXProblem", "run_simulations", "generate_equinox"]
