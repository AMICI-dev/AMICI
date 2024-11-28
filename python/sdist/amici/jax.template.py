import jax.numpy as jnp
from interpax import interp1d

from amici.jax.model import JAXModel


class JAXModel_TPL_MODEL_NAME(JAXModel):
    api_version = TPL_MODEL_API_VERSION

    def __init__(self):
        super().__init__()

    def _xdot(self, t, x, args):

        pk, tcl = args

        TPL_X_SYMS = x
        TPL_PK_SYMS = pk
        TPL_TCL_SYMS = tcl
        TPL_W_SYMS = self._w(t, x, pk, tcl)

        TPL_XDOT_EQ

        return TPL_XDOT_RET

    def _w(self, t, x, pk, tcl):

        TPL_X_SYMS = x
        TPL_PK_SYMS = pk
        TPL_TCL_SYMS = tcl

        TPL_W_EQ

        return TPL_W_RET

    def _x0(self, pk):

        TPL_PK_SYMS = pk

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

    def _tcl(self, x, pk):

        TPL_X_RDATA_SYMS = x
        TPL_PK_SYMS = pk

        TPL_TOTAL_CL_EQ

        return TPL_TOTAL_CL_RET

    def _y(self, t, x, pk, tcl):

        TPL_X_SYMS = x
        TPL_PK_SYMS = pk
        TPL_W_SYMS = self._w(t, x, pk, tcl)

        TPL_Y_EQ

        return TPL_Y_RET

    def _sigmay(self, y, pk):
        TPL_PK_SYMS = pk

        TPL_Y_SYMS = y

        TPL_SIGMAY_EQ

        return TPL_SIGMAY_RET


    def _nllh(self, t, x, pk, tcl, my, iy):
        y = self._y(t, x, pk, tcl)
        TPL_Y_SYMS = y
        TPL_SIGMAY_SYMS = self._sigmay(y, pk)

        TPL_JY_EQ

        return TPL_JY_RET.at[iy].get()

    @property
    def observable_ids(self):
        return TPL_Y_IDS

    @property
    def state_ids(self):
        return TPL_X_IDS

    @property
    def parameter_ids(self):
        return TPL_PK_IDS
