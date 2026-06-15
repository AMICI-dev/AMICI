"""Robertson stiff chemistry.

  dy1/dt = -k1*y1 + k2*y2*y3
  dy2/dt =  k1*y1 - k2*y2*y3 - k3*y2^2
  dy3/dt =  k3*y2^2

Exercises the adaptive step-size controller under stiffness.
"""

from pathlib import Path

import jax.numpy as jnp
from amici.sim.jax.model import JAXModel, safe_div, safe_log

DEFAULT_P = jnp.array([0.04, 1e4, 3e7])  # k1, k2, k3 (canonical values)
TS_DYN = jnp.array([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0])
# y1 starts near 1 and decays slowly — measurements slightly off
MY = jnp.array([1.0, 1.0, 1.0, 0.99, 0.97, 0.91, 0.71, 0.45])
IYS = jnp.zeros(len(TS_DYN), dtype=jnp.int32)
IY_TRAFOS = jnp.zeros(len(TS_DYN), dtype=jnp.int32)
OPS = jnp.zeros((len(TS_DYN), 0))
NPS = jnp.zeros((len(TS_DYN), 0))
SIGMA = 0.05


class Robertson(JAXModel):
    """Stiff Robertson chemistry — stresses diffrax adaptive step control."""

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
        y1, y2, y3 = x
        k1, k2, k3 = p
        dy1 = -k1 * y1 + k2 * y2 * y3
        dy2 = k1 * y1 - k2 * y2 * y3 - k3 * y2**2
        dy3 = k3 * y2**2
        return jnp.array([dy1, dy2, dy3])

    def _w(self, t, x, p, tcl, h):
        return jnp.array([])

    def _x0(self, t, p):
        return jnp.array([1.0, 0.0, 0.0])

    def _x_solver(self, x):
        return x

    def _x_rdata(self, x, tcl):
        return x

    def _tcl(self, x, p):
        return jnp.array([])

    def _y(self, t, x, p, tcl, h, op):
        y1, _y2, _y3 = x
        return jnp.array([y1])

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
        return ("y1", "y2", "y3")

    @property
    def observable_ids(self):
        return ("y1_obs",)

    @property
    def parameter_ids(self):
        return ("k1", "k2", "k3")

    @property
    def expression_ids(self):
        return ()
