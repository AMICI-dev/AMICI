"""Simple two-state system with a unique asymptotic steady state.

  dA/dt = k_prod - k_deg_A * A       →  A* = k_prod / k_deg_A   = 1.0
  dB/dt = k_deg_A * A - k_deg_B * B  →  B* = k_prod / k_deg_B   = 2.0

Used for two test operations:

* preeq  — call preequilibrate_condition() to reach (A*=1, B*=2) from (0,0).
           Stresses the eq() function.

* posteq — call simulate_condition_unjitted() with ts_posteq=[inf],
           so eq() is called after the dynamic simulation.
"""

from pathlib import Path

import jax.numpy as jnp
from amici.sim.jax.model import JAXModel, safe_div, safe_log

DEFAULT_P = jnp.array([1.0, 1.0, 0.5])  # k_prod, k_deg_A, k_deg_B

# Post-eq configuration: simulate for a short time then equilibrate
TS_DYN_POSTEQ = jnp.array([0.5, 1.0, 2.0, 3.0])
TS_POSTEQ = jnp.array([jnp.inf])  # one observable measured at steady state
_NT_POSTEQ = len(TS_DYN_POSTEQ) + len(TS_POSTEQ)
MY_POSTEQ = jnp.concatenate(
    [jnp.array([0.4, 0.6, 0.8, 0.9]), jnp.array([2.0])]
)  # B measurements (observable = B); last entry at post-eq steady state
IYS_POSTEQ = jnp.zeros(_NT_POSTEQ, dtype=jnp.int32)
IY_TRAFOS_POSTEQ = jnp.zeros(_NT_POSTEQ, dtype=jnp.int32)
OPS_POSTEQ = jnp.zeros((_NT_POSTEQ, 0))
NPS_POSTEQ = jnp.zeros((_NT_POSTEQ, 0))
SIGMA = 0.1


class Equilibration(JAXModel):
    """Two-state production/degradation system used for preeq and posteq tests."""

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
        A, B = x
        k_prod, k_deg_A, k_deg_B = p
        dA = k_prod - k_deg_A * A
        dB = k_deg_A * A - k_deg_B * B
        return jnp.array([dA, dB])

    def _w(self, t, x, p, tcl, h):
        return jnp.array([])

    def _x0(self, t, p):
        return jnp.array([0.0, 0.0])  # start far from steady state

    def _x_solver(self, x):
        return x

    def _x_rdata(self, x, tcl):
        return x

    def _tcl(self, x, p):
        return jnp.array([])

    def _y(self, t, x, p, tcl, h, op):
        _A, B = x
        return jnp.array([B])  # observable = B

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
        return ("B_obs",)

    @property
    def parameter_ids(self):
        return ("k_prod", "k_deg_A", "k_deg_B")

    @property
    def expression_ids(self):
        return ()
