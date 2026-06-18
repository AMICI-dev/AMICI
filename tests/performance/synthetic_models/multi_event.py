"""Three sequential impulses with a constant-zero ODE.

  dx/dt = 0
  At t=2: x += delta1
  At t=4: x += delta2
  At t=6: x += delta3

Uses three event pairs (6 heaviside variables) to stress the while_loop in
solve() and _check_cascading_events.
"""

from pathlib import Path

import jax.numpy as jnp
from amici.sim.jax.model import JAXModel, safe_div, safe_log

DEFAULT_P = jnp.array([0.5, -0.3, 0.8])  # delta1, delta2, delta3
TS_DYN = jnp.array([1.0, 2.5, 3.5, 4.5, 5.5, 6.5, 8.0])
# x=0, after t=2: x=0.5, after t=4: x=0.2, after t=6: x=1.0
MY = jnp.array([0.0, 0.5, 0.5, 0.2, 0.2, 1.0, 1.0])
IYS = jnp.zeros(len(TS_DYN), dtype=jnp.int32)
IY_TRAFOS = jnp.zeros(len(TS_DYN), dtype=jnp.int32)
OPS = jnp.zeros((len(TS_DYN), 0))
NPS = jnp.zeros((len(TS_DYN), 0))
SIGMA = 0.05


class MultiEvent(JAXModel):
    """Three sequential impulses — stresses the while_loop in solve()."""

    api_version = "0.0.4"

    def __init__(self):
        self.jax_py_file = Path(__file__).resolve()
        self.nns = {}
        self._array_inputs = {}
        self._array_input_index = jnp.int32(0)
        self.parameters = DEFAULT_P
        super().__init__()

    def _xdot(self, t, x, args):
        return jnp.array([0.0])

    def _w(self, t, x, p, tcl, h):
        return jnp.array([])

    def _x0(self, t, p):
        return jnp.array([0.0])

    def _x_solver(self, x):
        return x

    def _x_rdata(self, x, tcl):
        return x

    def _tcl(self, x, p):
        return jnp.array([])

    def _y(self, t, x, p, tcl, h, op):
        (state,) = x
        return jnp.array([state])

    def _sigmay(self, y, p, np):
        return jnp.array([SIGMA])

    def _nllh(self, t, x, p, tcl, h, my, iy, op, np):
        y = self._y(t, x, p, tcl, h, op)
        sigma = self._sigmay(y, p, np).at[iy].get()
        y_val = y.at[iy].get()
        return 0.5 * safe_log(2 * jnp.pi * sigma**2) + safe_div(
            0.5 * (y_val - my) ** 2, sigma**2
        )

    # ── Events ─────────────────────────────────────────────────────────────
    # Three pairs at t=2, t=4, t=6.  Only the even-indexed events fire.

    def _known_discs(self, p):
        return jnp.array([2.0, 4.0, 6.0])

    def _root_cond_fn(self, t, y, args, **_):
        return jnp.array(
            [
                t - 2.0,
                2.0 - t,  # pair 0: fires at t=2
                t - 4.0,
                4.0 - t,  # pair 1: fires at t=4
                t - 6.0,
                6.0 - t,  # pair 2: fires at t=6
            ]
        )

    def _delta_x(self, y, p, tcl):
        # delta_x shape: (n_events * n_solver_states,) = (6*1,) = (6,)
        # Odd-indexed events (1, 3, 5) never fire.
        d1, d2, d3 = p
        return jnp.array([d1, 0.0, d2, 0.0, d3, 0.0])

    @property
    def event_initial_values(self):
        return jnp.array([True, True, True, True, True, True])

    @property
    def n_events(self):
        return 6

    @property
    def state_ids(self):
        return ("x",)

    @property
    def observable_ids(self):
        return ("x_obs",)

    @property
    def parameter_ids(self):
        return ("delta1", "delta2", "delta3")

    @property
    def expression_ids(self):
        return ()
