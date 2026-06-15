"""Linear decay: dx/dt = -k*x.  Analytical solution x(t) = x0*exp(-k*t)."""

from pathlib import Path

import jax.numpy as jnp
from amici.sim.jax.model import JAXModel, safe_div, safe_log

# Fixed test configuration ──────────────────────────────────────────────────
DEFAULT_P = jnp.array([0.1, 1.0])  # k=0.1, x0=1.0
TS_DYN = jnp.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
# Measurements slightly off the true trajectory so llh is finite and non-trivial
MY = jnp.array([0.91, 0.83, 0.75, 0.68, 0.62, 0.56, 0.51, 0.46, 0.42, 0.38])
IYS = jnp.zeros(len(TS_DYN), dtype=jnp.int32)
IY_TRAFOS = jnp.zeros(len(TS_DYN), dtype=jnp.int32)
OPS = jnp.zeros((len(TS_DYN), 0))
NPS = jnp.zeros((len(TS_DYN), 0))
SIGMA = 0.1


class LinearDecay(JAXModel):
    """Canary model — any regression here indicates a fundamental problem."""

    api_version = "0.0.4"

    def __init__(self):
        self.jax_py_file = Path(__file__).resolve()
        self.nns = {}
        self._array_inputs = {}
        self._array_input_index = jnp.int32(0)
        self.parameters = DEFAULT_P
        super().__init__()

    # ── ODE ────────────────────────────────────────────────────────────────

    def _xdot(self, t, x, args):
        p, tcl, h = args
        (state,) = x
        k, _x0 = p
        return jnp.array([-k * state])

    def _w(self, t, x, p, tcl, h):
        return jnp.array([])

    # ── State / initial conditions ─────────────────────────────────────────

    def _x0(self, t, p):
        _k, x0 = p
        return jnp.array([x0])

    def _x_solver(self, x):
        return x

    def _x_rdata(self, x, tcl):
        return x

    def _tcl(self, x, p):
        return jnp.array([])

    # ── Observables / likelihood ───────────────────────────────────────────

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

    # ── Events (none) ──────────────────────────────────────────────────────

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

    # ── Identifiers ───────────────────────────────────────────────────────

    @property
    def state_ids(self):
        return ("x",)

    @property
    def observable_ids(self):
        return ("x_obs",)

    @property
    def parameter_ids(self):
        return ("k", "x0")

    @property
    def expression_ids(self):
        return ()
