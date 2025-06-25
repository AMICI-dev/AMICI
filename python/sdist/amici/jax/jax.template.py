# ruff: noqa: F401, F821, F841
import jax.numpy as jnp
import jaxtyping as jt
from interpax import interp1d
from pathlib import Path
from jax.numpy import inf as oo
from jax.numpy import nan as nan

from amici.jax.model import JAXModel, safe_log, safe_div


class JAXModel_TPL_MODEL_NAME(JAXModel):
    api_version = TPL_MODEL_API_VERSION

    def __init__(self):
        self.jax_py_file = Path(__file__).resolve()
        self.parameters = TPL_P_VALUES
        super().__init__()

    def _xdot(self, t, x, args):
        p, tcl = args

        TPL_X_SYMS = x
        TPL_P_SYMS = p
        TPL_TCL_SYMS = tcl
        TPL_W_SYMS = self._w(t, x, p, tcl)

        TPL_XDOT_EQ

        return TPL_XDOT_RET

    def _w(self, t, x, p, tcl):
        TPL_X_SYMS = x
        TPL_P_SYMS = p
        TPL_TCL_SYMS = tcl

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

    def _y(self, t, x, p, tcl, op):
        TPL_X_SYMS = x
        TPL_P_SYMS = p
        TPL_W_SYMS = self._w(t, x, p, tcl)
        TPL_OP_SYMS = op

        TPL_Y_EQ

        return TPL_Y_RET

    def _sigmay(self, y, p, np):
        TPL_P_SYMS = p

        TPL_Y_SYMS = y
        TPL_NP_SYMS = np

        TPL_SIGMAY_EQ

        return TPL_SIGMAY_RET

    def _nllh(self, t, x, p, tcl, my, iy, op, np):
        y = self._y(t, x, p, tcl, op)
        if not y.size:
            return jnp.array(0.0)

        TPL_Y_SYMS = y
        TPL_SIGMAY_SYMS = self._sigmay(y, p, np)

        TPL_JY_EQ

        return TPL_JY_RET.at[iy].get()

    def _known_discs(self, p):
        TPL_P_SYMS = p

        return TPL_ROOTS

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
