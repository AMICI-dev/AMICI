# ruff: noqa: F401, F821, F841
from pathlib import Path

import equinox as eqx
import jax.numpy as jnp
import jax.random as jr
import jaxtyping as jt
from interpax import interp1d
from jax.numpy import inf as oo
from jax.numpy import nan as nan

from amici import _module_from_path
from amici.jax.model import JAXModel, safe_div, safe_log

TPL_NET_IMPORTS


class JAXModel_TPL_MODEL_NAME(JAXModel):
    api_version = TPL_MODEL_API_VERSION

    def __init__(self):
        self.jax_py_file = Path(__file__).resolve()
        self.nns = {TPL_NETS}
        self.parameters = TPL_P_VALUES
        super().__init__()

    def _xdot(self, t, x, args):
        p, tcl, h = args

        TPL_X_SYMS = x
        TPL_P_SYMS = p
        TPL_TCL_SYMS = tcl
        TPL_IH_SYMS = h
        TPL_W_SYMS = self._w(t, x, p, tcl, h)

        TPL_XDOT_EQ

        return TPL_XDOT_RET

    def _w(self, t, x, p, tcl, h):
        TPL_X_SYMS = x
        TPL_P_SYMS = p
        TPL_TCL_SYMS = tcl
        TPL_IH_SYMS = h

        TPL_W_EQ

        return TPL_W_RET

    def _x0(self, t, p):
        TPL_P_SYMS = p

        TPL_X0_EQ

        return TPL_X0_RET

    def _x_solver(self, x):
        TPL_X_RDATA_SYMS = x

        TPL_X_SOLVER_EQ

        return TPL_X_SOLVER_RET

    def _x_rdata(self, x, tcl):
        TPL_X_SYMS = x
        TPL_TCL_SYMS = tcl

        TPL_X_RDATA_EQ

        return TPL_X_RDATA_RET

    def _tcl(self, x, p):
        TPL_X_RDATA_SYMS = x
        TPL_P_SYMS = p

        TPL_TOTAL_CL_EQ

        return TPL_TOTAL_CL_RET

    def _y(self, t, x, p, tcl, h, op):
        TPL_X_SYMS = x
        TPL_P_SYMS = p
        TPL_W_SYMS = self._w(t, x, p, tcl, h)
        TPL_OP_SYMS = op

        TPL_Y_EQ

        return TPL_Y_RET

    def _sigmay(self, y, p, np):
        TPL_P_SYMS = p

        TPL_Y_SYMS = y
        TPL_NP_SYMS = np

        TPL_SIGMAY_EQ

        return TPL_SIGMAY_RET

    def _nllh(self, t, x, p, tcl, h, my, iy, op, np):
        y = self._y(t, x, p, tcl, h, op)
        if not y.size:
            return jnp.array(0.0)

        TPL_Y_SYMS = y
        TPL_SIGMAY_SYMS = self._sigmay(y, p, np)

        TPL_JY_EQ

        return TPL_JY_RET.at[iy].get()

    def _known_discs(self, p):
        TPL_P_SYMS = p

        return TPL_ROOTS

    def _root_cond_fn(self, t, y, args, **_):
        p, tcl, h = args

        TPL_X_SYMS = y
        TPL_P_SYMS = p
        TPL_TCL_SYMS = tcl
        TPL_IH_SYMS = h
        TPL_W_SYMS = self._w(t, y, p, tcl, h)

        TPL_IROOT_EQ
        TPL_EROOT_EQ

        return jnp.hstack((TPL_IROOT_RET, TPL_EROOT_RET))
    
    def _delta_x(self, y, p, tcl):
        TPL_X_SYMS = y
        TPL_P_SYMS = p
        TPL_TCL_SYMS = tcl

        TPL_X_OLD_EQ
        
        TPL_DELTAX_EQ

        return TPL_DELTAX_RET
    
    @property
    def event_initial_values(self):
        return TPL_EVENT_INITIAL_VALUES

    @property
    def n_events(self):
        return TPL_N_IEVENTS + TPL_N_EEVENTS

    @property
    def observable_ids(self):
        return TPL_Y_IDS

    @property
    def state_ids(self):
        return TPL_X_RDATA_IDS

    @property
    def parameter_ids(self):
        return TPL_P_IDS

    @property
    def expression_ids(self):
        return TPL_W_IDS


Model = JAXModel_TPL_MODEL_NAME
