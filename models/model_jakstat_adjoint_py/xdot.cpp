#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"
#include "xdot.h"

namespace amici {
namespace model_model_jakstat_adjoint_py {

void xdot_model_jakstat_adjoint_py(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dSTATdt = (-Omega_cyt*STAT*p1*u + Omega_nuc*nSTAT5*p4)/Omega_cyt;  // xdot[0]
    dpSTATdt = STAT*p1*u - 2*p2*std::pow(pSTAT, 2);  // xdot[1]
    dpSTAT_pSTATdt = p2*std::pow(pSTAT, 2) - p3*pSTAT_pSTAT;  // xdot[2]
    dnpSTAT_npSTATdt = (Omega_cyt*p3*pSTAT_pSTAT - Omega_nuc*npSTAT_npSTAT*p4)/Omega_nuc;  // xdot[3]
    dnSTAT1dt = -p4*(nSTAT1 - 2*npSTAT_npSTAT);  // xdot[4]
    dnSTAT2dt = p4*(nSTAT1 - nSTAT2);  // xdot[5]
    dnSTAT3dt = p4*(nSTAT2 - nSTAT3);  // xdot[6]
    dnSTAT4dt = p4*(nSTAT3 - nSTAT4);  // xdot[7]
    dnSTAT5dt = p4*(nSTAT4 - nSTAT5);  // xdot[8]
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
