"""Reversible conversion A ⇌ B with one conservation law.

  Full states: [A, B],  Conservation: A + B = total (= 1)
  Solver state: [A] only.  B reconstructed as total - A.

  dA/dt = -k1*A + k2*B  =  -(k1+k2)*A + k2*total

Exercises _tcl(), _x_solver(), _x_rdata().
"""

from pathlib import Path

import jax.numpy as jnp
from amici.sim.jax.model import JAXModel, safe_div, safe_log

DEFAULT_P = jnp.array([0.5, 0.2])  # k1 (A->B), k2 (B->A)
# Steady state: A* = k2/(k1+k2) = 0.2/0.7 ≈ 0.286, B* ≈ 0.714
TS_DYN = jnp.array([0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0])
MY = jnp.array([0.75, 0.60, 0.45, 0.38, 0.33, 0.30, 0.29, 0.29, 0.29, 0.29])
IYS = jnp.zeros(len(TS_DYN), dtype=jnp.int32)
IY_TRAFOS = jnp.zeros(len(TS_DYN), dtype=jnp.int32)
OPS = jnp.zeros((len(TS_DYN), 0))
NPS = jnp.zeros((len(TS_DYN), 0))
SIGMA = 0.05


class ConservationLaw(JAXModel):
    """Two-state reversible system with one conservation law."""

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
        (A,) = x  # solver state contains only A
        (total,) = tcl
        k1, k2 = p
        B = total - A
        dA = -k1 * A + k2 * B
        return jnp.array([dA])

    def _w(self, t, x, p, tcl, h):
        return jnp.array([])

    def _x0(self, t, p):
        return jnp.array([1.0, 0.0])  # full state [A, B]

    def _x_solver(self, x):
        A, _B = x
        return jnp.array([A])

    def _x_rdata(self, x, tcl):
        (A,) = x
        (total,) = tcl
        return jnp.array([A, total - A])

    def _tcl(self, x, p):
        A, B = x
        return jnp.array([A + B])

    def _y(self, t, x, p, tcl, h, op):
        (A,) = x
        return jnp.array([A])

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
        return ("A", "B")

    @property
    def observable_ids(self):
        return ("A_obs",)

    @property
    def parameter_ids(self):
        return ("k1", "k2")

    @property
    def expression_ids(self):
        return ()
