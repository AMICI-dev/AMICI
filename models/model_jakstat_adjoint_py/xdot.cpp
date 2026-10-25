#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_jakstat_adjoint_py {

void xdot_model_jakstat_adjoint_py(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    const realtype STAT_ = x[0];
    const realtype pSTAT_ = x[1];
    const realtype pSTAT_pSTAT_ = x[2];
    const realtype npSTAT_npSTAT_ = x[3];
    const realtype nSTAT1_ = x[4];
    const realtype nSTAT2_ = x[5];
    const realtype nSTAT3_ = x[6];
    const realtype nSTAT4_ = x[7];
    const realtype nSTAT5_ = x[8];
    const realtype p1_ = p[0];
    const realtype p2_ = p[1];
    const realtype p3_ = p[2];
    const realtype p4_ = p[3];
    const realtype Omega_cyt_ = k[0];
    const realtype Omega_nuc_ = k[1];
    const realtype u_ = w[0];

    realtype &dSTATdt_ = xdot[0];
    realtype &dpSTATdt_ = xdot[1];
    realtype &dpSTAT_pSTATdt_ = xdot[2];
    realtype &dnpSTAT_npSTATdt_ = xdot[3];
    realtype &dnSTAT1dt_ = xdot[4];
    realtype &dnSTAT2dt_ = xdot[5];
    realtype &dnSTAT3dt_ = xdot[6];
    realtype &dnSTAT4dt_ = xdot[7];
    realtype &dnSTAT5dt_ = xdot[8];
    dSTATdt_ = (-Omega_cyt_*STAT_*p1_*u_ + Omega_nuc_*nSTAT5_*p4_)/Omega_cyt_;  // xdot[0]
    dpSTATdt_ = STAT_*p1_*u_ - 2*p2_*std::pow(pSTAT_, 2);  // xdot[1]
    dpSTAT_pSTATdt_ = p2_*std::pow(pSTAT_, 2) - p3_*pSTAT_pSTAT_;  // xdot[2]
    dnpSTAT_npSTATdt_ = (Omega_cyt_*p3_*pSTAT_pSTAT_ - Omega_nuc_*npSTAT_npSTAT_*p4_)/Omega_nuc_;  // xdot[3]
    dnSTAT1dt_ = -p4_*(nSTAT1_ - 2*npSTAT_npSTAT_);  // xdot[4]
    dnSTAT2dt_ = p4_*(nSTAT1_ - nSTAT2_);  // xdot[5]
    dnSTAT3dt_ = p4_*(nSTAT2_ - nSTAT3_);  // xdot[6]
    dnSTAT4dt_ = p4_*(nSTAT3_ - nSTAT4_);  // xdot[7]
    dnSTAT5dt_ = p4_*(nSTAT4_ - nSTAT5_);  // xdot[8]
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
