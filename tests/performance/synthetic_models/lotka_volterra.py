"""Lotka-Volterra predator-prey oscillator.

  dx/dt = alpha*x - beta*x*y   (prey)
  dy/dt = -gamma*y + delta*x*y  (predator)

Exercises long-horizon non-stiff integration.
"""

from pathlib import Path

import jax.numpy as jnp
from amici.sim.jax.model import JAXModel, safe_div, safe_log

DEFAULT_P = jnp.array([1.5, 1.0, 3.0, 1.0])  # alpha, beta, gamma, delta
TS_DYN = jnp.linspace(0.5, 12.0, 15)
# Prey starts at 10, oscillates — use approximate reference values
MY = jnp.array(
    [
        9.0,
        4.0,
        2.0,
        4.0,
        10.0,
        15.0,
        10.0,
        4.0,
        2.0,
        4.0,
        10.0,
        15.0,
        10.0,
        4.0,
        2.0,
    ]
)
IYS = jnp.zeros(len(TS_DYN), dtype=jnp.int32)
IY_TRAFOS = jnp.zeros(len(TS_DYN), dtype=jnp.int32)
OPS = jnp.zeros((len(TS_DYN), 0))
NPS = jnp.zeros((len(TS_DYN), 0))
SIGMA = 2.0


class LotkaVolterra(JAXModel):
    """Non-stiff oscillator — exercises long-horizon integration."""

    api_version = "0.0.4"

    def __init__(self):
        self.jax_py_file = Path(__file__).resolve()
        self.nns = {}
        self._array_inputs = {}
        self._array_input_index = jnp.int32(0)
        self.parameters = DEFAULT_P
        super().__init__()

    def _xdot(self, t, x, args):
        p, tcl, h = args
        prey, pred = x
        alpha, beta, gamma, delta = p
        dprey = alpha * prey - beta * prey * pred
        dpred = -gamma * pred + delta * prey * pred
        return jnp.array([dprey, dpred])

    def _w(self, t, x, p, tcl, h):
        return jnp.array([])

    def _x0(self, t, p):
        return jnp.array([10.0, 5.0])

    def _x_solver(self, x):
        return x

    def _x_rdata(self, x, tcl):
        return x

    def _tcl(self, x, p):
        return jnp.array([])

    def _y(self, t, x, p, tcl, h, op):
        prey, _pred = x
        return jnp.array([prey])

    def _sigmay(self, y, p, np):
        return jnp.array([SIGMA])

    def _nllh(self, t, x, p, tcl, h, my, iy, op, np):
        y = self._y(t, x, p, tcl, h, op)
        sigma = self._sigmay(y, p, np).at[iy].get()
        y_val = y.at[iy].get()
        return 0.5 * safe_log(2 * jnp.pi * sigma**2) + safe_div(
            0.5 * (y_val - my) ** 2, sigma**2
        )

    def _known_discs(self, p):
        return jnp.array([])

    def _root_cond_fn(self, t, y, args, **_):
        return jnp.array([])

    def _delta_x(self, y, p, tcl):
        return jnp.array([])

    @property
    def event_initial_values(self):
        return jnp.array([])

    @property
    def n_events(self):
        return 0

    @property
    def state_ids(self):
        return ("prey", "predator")

    @property
    def observable_ids(self):
        return ("prey_obs",)

    @property
    def parameter_ids(self):
        return ("alpha", "beta", "gamma", "delta")

    @property
    def expression_ids(self):
        return ()
