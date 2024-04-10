import jax.numpy as jnp

from amici.jax import JAXModel


class JAXModel_TPL_MODEL_NAME(JAXModel):
    def __init__(self):
        super().__init__()

    def xdot(self, t, x, args):

        p, k, tcl = args

        TPL_X_SYMS = x
        TPL_P_SYMS = p
        TPL_K_SYMS = k
        TPL_TCL_SYMS = tcl
        TPL_W_SYMS = self._w(t, x, p, k, tcl)

TPL_XDOT_EQ

        return TPL_XDOT_RET

    def _w(self, t, x, p, k, tcl):

        TPL_X_SYMS = x
        TPL_P_SYMS = p
        TPL_K_SYMS = k
        TPL_TCL_SYMS = tcl

TPL_W_EQ

        return TPL_W_RET

    def x0(self, p, k):

        TPL_P_SYMS = p
        TPL_K_SYMS = k

TPL_X0_EQ

        return TPL_X0_RET

    def x_solver(self, x):

        TPL_X_RDATA_SYMS = x

TPL_X_SOLVER_EQ

        return TPL_X_SOLVER_RET

    def x_rdata(self, x, tcl):

        TPL_X_SYMS = x
        TPL_TCL_SYMS = tcl

TPL_X_RDATA_EQ

        return TPL_X_RDATA_RET

    def tcl(self, x, p, k):

        TPL_X_RDATA_SYMS = x
        TPL_P_SYMS = p
        TPL_K_SYMS = k

TPL_TOTAL_CL_EQ

        return TPL_TOTAL_CL_RET

    def y(self, t, x, p, k, tcl):

        TPL_X_SYMS = x
        TPL_P_SYMS = p
        TPL_K_SYMS = k
        TPL_W_SYMS = self._w(t, x, p, k, tcl)

TPL_Y_EQ

        return TPL_Y_RET

    def sigmay(self, y, p, k):
        TPL_Y_SYMS = y
        TPL_P_SYMS = p
        TPL_K_SYMS = k

TPL_SIGMAY_EQ

        return TPL_SIGMAY_RET

    def Jy(self, y, my, sigmay):
        TPL_Y_SYMS = y
        TPL_MY_SYMS = my
        TPL_SIGMAY_SYMS = sigmay

TPL_JY_EQ

        return TPL_JY_RET
