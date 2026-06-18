"""Single-event model: constant growth with one impulse.

  dx/dt = r
  At t = t_event: x += delta

Uses one event pair (2 heaviside variables) to test event detection and
the _run_segment → event → _run_segment code path in solve().
"""

from pathlib import Path

import jax.numpy as jnp
from amici.sim.jax.model import JAXModel, safe_div, safe_log

# k=r, p[1]=delta, p[2]=t_event
DEFAULT_P = jnp.array([0.1, 0.5, 5.0])  # r=0.1, delta=0.5, t_event=5.0
# Time points spanning both sides of the event at t=5
TS_DYN = jnp.array([1.0, 3.0, 4.9, 5.1, 7.0, 10.0])
# x starts at 0: x(t<5) = r*t, x(t>=5) = r*t + delta
MY = jnp.array([0.1, 0.3, 0.49, 1.01, 1.2, 1.5])
IYS = jnp.zeros(len(TS_DYN), dtype=jnp.int32)
IY_TRAFOS = jnp.zeros(len(TS_DYN), dtype=jnp.int32)
OPS = jnp.zeros((len(TS_DYN), 0))
NPS = jnp.zeros((len(TS_DYN), 0))
SIGMA = 0.05


class SingleEvent(JAXModel):
    """Constant growth with one impulse — tests event detection / _handle_event."""

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
        r, _delta, _t_event = p
        return jnp.array([r])

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
    # One event pair: event 0 fires at t_event (rising root), event 1 never fires.

    def _known_discs(self, p):
        _r, _delta, t_event = p
        return jnp.array([t_event])

    def _root_cond_fn(self, t, y, args, **_):
        p, tcl, h = args
        _r, _delta, t_event = p
        eroot0 = t - t_event  # rises through 0 at t_event  → fires
        eroot1 = (
            t_event - t
        )  # falls through 0 at t_event  → never fires (h[1]=1 from init)
        return jnp.array([eroot0, eroot1])

    def _delta_x(self, y, p, tcl):
        # delta_x shape: (n_events * n_solver_states,) = (2*1,) = (2,)
        # Event 0 fires: x += delta
        # Event 1 never fires: 0
        _r, delta, _t_event = p
        return jnp.array([delta, 0.0])

    @property
    def event_initial_values(self):
        # h0 = [True, True] → after _initialise_heaviside_variables → h = [0, 1]
        return jnp.array([True, True])

    @property
    def n_events(self):
        return 2

    @property
    def state_ids(self):
        return ("x",)

    @property
    def observable_ids(self):
        return ("x_obs",)

    @property
    def parameter_ids(self):
        return ("r", "delta", "t_event")

    @property
    def expression_ids(self):
        return ()
