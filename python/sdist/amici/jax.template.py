import jax.numpy as jnp
from interpax import interp1d

from amici.jax import JAXModel


class JAXModel_TPL_MODEL_NAME(JAXModel):
    def __init__(self):
        super().__init__()

    @staticmethod
    def xdot(t, x, args):

        p, k, tcl = args

        TPL_X_SYMS = x
        TPL_P_SYMS = p
        TPL_K_SYMS = k
        TPL_TCL_SYMS = tcl
        TPL_W_SYMS = JAXModel_TPL_MODEL_NAME._w(t, x, p, k, tcl)

TPL_XDOT_EQ

        return TPL_XDOT_RET

    @staticmethod
    def _w(t, x, p, k, tcl):

        TPL_X_SYMS = x
        TPL_P_SYMS = p
        TPL_K_SYMS = k
        TPL_TCL_SYMS = tcl

TPL_W_EQ

        return TPL_W_RET

    @staticmethod
    def x0(p, k):

        TPL_P_SYMS = p
        TPL_K_SYMS = k

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
    def tcl(x, p, k):

        TPL_X_RDATA_SYMS = x
        TPL_P_SYMS = p
        TPL_K_SYMS = k

TPL_TOTAL_CL_EQ

        return TPL_TOTAL_CL_RET

    @staticmethod
    def y(t, x, p, k, tcl):

        TPL_X_SYMS = x
        TPL_P_SYMS = p
        TPL_K_SYMS = k
        TPL_W_SYMS = JAXModel_TPL_MODEL_NAME._w(t, x, p, k, tcl)

TPL_Y_EQ

        return TPL_Y_RET

    @staticmethod
    def sigmay(y, p, k):
        TPL_Y_SYMS = y
        TPL_P_SYMS = p
        TPL_K_SYMS = k

TPL_SIGMAY_EQ

        return TPL_SIGMAY_RET

    @staticmethod
    def Jy(y, my, sigmay):
        TPL_Y_SYMS = y
        TPL_MY_SYMS = my
        TPL_SIGMAY_SYMS = sigmay

TPL_JY_EQ

        return TPL_JY_RET

    @property
    def parameter_ids(self):
        return TPL_P_IDS

    @property
    def fixed_parameter_ids(self):
        return TPL_K_IDS

    @property
    def observable_ids(self):
        return TPL_Y_IDS

    @property
    def state_ids(self):
        return TPL_X_IDS
