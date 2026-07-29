# ruff: noqa: F401, F821, F841
from pathlib import Path

import equinox as eqx  # noqa: F401
import jax.numpy as jnp
import jax.random as jr  # noqa: F401
import jaxtyping as jt  # noqa: F401
from jax.numpy import inf as oo  # noqa: F401
from jax.numpy import nan as nan  # noqa: F401

from amici import _module_from_path  # noqa: F401
from amici.sim.jax.model import JAXModel, safe_div, safe_log  # noqa: F401




class JAXModel_test_arbitrary_placeholder_jax(JAXModel):
    api_version = '0.0.4'

    def __init__(self):
        self.jax_py_file = Path(__file__).resolve()
        self.nns = {}
        self.parameters = jnp.array([])
        self._array_inputs = {}
        self._array_input_index = jnp.int32(0)
        super().__init__()

    def _xdot(self, t, x, args):
        p, tcl, h = args

        xa, xb,  = x
        _ = p
        _ = tcl
        _ = h
        _ = self._w(t, x, p, tcl, h)

        dxadt = safe_div(-xa, 10)
        dxbdt = safe_div(xb, 20)

        return jnp.array([dxadt, dxbdt])

    def _w(self, t, x, p, tcl, h):
        xa, xb,  = x
        _ = p
        _ = tcl
        _ = h

        

        return jnp.array([])

    def _x0(self, t, p):
        _ = p

        x00 = 1.00000000000000
        x01 = 2.00000000000000

        return jnp.array([x00, x01])

    def _x_solver(self, x):
        xa, xb,  = x

        x_solver0 = xa
        x_solver1 = xb

        return jnp.array([x_solver0, x_solver1])

    def _x_rdata(self, x, tcl):
        xa, xb,  = x
        _ = tcl

        xa = xa
        xb = xb

        return jnp.array([xa, xb])

    def _tcl(self, x, p):
        xa, xb,  = x
        _ = p

        

        return jnp.array([])

    def _y(self, t, x, p, tcl, h, op):
        xa, xb,  = x
        _ = p
        _ = self._w(t, x, p, tcl, h)
        _ = op

        obs_a = xa
        obs_b = xb

        return jnp.array([obs_a, obs_b])

    def _sigmay(self, y, p, np):
        _ = p

        obs_a, obs_b,  = y
        noiseParameter1,  = np

        sigma_obs_a = noiseParameter1*obs_a
        sigma_obs_b = noiseParameter1

        return jnp.array([sigma_obs_a, sigma_obs_b])

    def _nllh(self, t, x, p, tcl, h, my, iy, op, np):
        y = self._y(t, x, p, tcl, h, op)
        if not y.size:
            return jnp.array(0.0)

        obs_a, obs_b,  = y
        sigma_obs_a, sigma_obs_b,  = self._sigmay(y, p, np)

        _amici_cse_0 = sigma_obs_a**2
        _amici_cse_1 = 2*jnp.pi
        _amici_cse_2 = sigma_obs_b**2
        Jy0 = 0.5*safe_log(_amici_cse_0*_amici_cse_1) + safe_div(0.5*(-my + obs_a)**2, _amici_cse_0)
        Jy1 = 0.5*safe_log(_amici_cse_1*_amici_cse_2) + safe_div(0.5*(-my + obs_b)**2, _amici_cse_2)

        return jnp.array([Jy0, Jy1]).at[iy].get()

    def _known_discs(self, p):
        _ = p

        return jnp.array([])

    def _root_cond_fn(self, t, y, args, **_):
        p, tcl, h = args

        xa, xb,  = y
        _ = p
        _ = tcl
        _ = h
        _ = self._w(t, y, p, tcl, h)

        
        

        return jnp.hstack((jnp.array([]), jnp.array([])))

    def _delta_x(self, y, p, tcl):
        xa, xb,  = y
        _ = p
        _ = tcl
        # FIXME: workaround until state from event time is properly passed
        x_old0, x_old1,  = y

        

        return jnp.array([])

    @property
    def event_initial_values(self):
        return jnp.array([])

    @property
    def n_events(self):
        return 0 + 0

    @property
    def observable_ids(self):
        return "obs_a", "obs_b", 

    @property
    def state_ids(self):
        return "xa", "xb", 

    @property
    def parameter_ids(self):
        return tuple()

    @property
    def expression_ids(self):
        return tuple()


Model = JAXModel_test_arbitrary_placeholder_jax
