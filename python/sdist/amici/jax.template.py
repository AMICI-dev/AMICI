import jax.numpy as jnp
from interpax import interp1d

from amici.jax.model import JAXModel


class JAXModel_TPL_MODEL_NAME(JAXModel):
    def __init__(self):
        super().__init__()

    @staticmethod
    def xdot(t, x, args):

        pk, tcl = args

        TPL_X_SYMS = x
        TPL_PK_SYMS = pk
        TPL_TCL_SYMS = tcl
        TPL_W_SYMS = JAXModel_TPL_MODEL_NAME._w(t, x, pk, tcl)

TPL_XDOT_EQ

        return TPL_XDOT_RET

    @staticmethod
    def _w(t, x, pk, tcl):

        TPL_X_SYMS = x
        TPL_PK_SYMS = pk
        TPL_TCL_SYMS = tcl

TPL_W_EQ

        return TPL_W_RET

    @staticmethod
    def x0(pk):

        TPL_PK_SYMS = pk

TPL_X0_EQ

        return TPL_X0_RET

    @staticmethod
    def x_solver(x):

        TPL_X_RDATA_SYMS = x

TPL_X_SOLVER_EQ

        return TPL_X_SOLVER_RET

    @staticmethod
    def x_rdata(x, tcl):

        TPL_X_SYMS = x
        TPL_TCL_SYMS = tcl

TPL_X_RDATA_EQ

        return TPL_X_RDATA_RET

    @staticmethod
    def tcl(x, pk):

        TPL_X_RDATA_SYMS = x
        TPL_PK_SYMS = pk

TPL_TOTAL_CL_EQ

        return TPL_TOTAL_CL_RET

    def y(self, t, x, pk, tcl):

        TPL_X_SYMS = x
        TPL_PK_SYMS = pk
        TPL_W_SYMS = self._w(t, x, pk, tcl)

TPL_Y_EQ

        return TPL_Y_RET

    def sigmay(self, y, pk):
        TPL_PK_SYMS = pk

        TPL_Y_SYMS = y

TPL_SIGMAY_EQ

        return TPL_SIGMAY_RET


    def llh(self, t, x, pk, tcl, my, iy):
        y = self.y(t, x, pk, tcl)
        TPL_Y_SYMS = y
        sigmay = self.sigmay(y, pk)
        TPL_SIGMAY_SYMS = sigmay

TPL_JY_EQ

        return jnp.array([
            TPL_JY_RET.at[iy].get(),
            y.at[iy].get(),
            sigmay.at[iy].get()
        ])

    @property
    def observable_ids(self):
        return TPL_Y_IDS

    @property
    def state_ids(self):
        return TPL_X_IDS

    @property
    def parameter_ids(self):
        return TPL_PK_IDS
