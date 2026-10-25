#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_jakstat_adjoint_py {

void dydx_model_jakstat_adjoint_py(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    const realtype init_STAT_ = p[4];
    const realtype scale_tSTAT_ = p[12];
    const realtype scale_pSTAT_ = p[13];

    dydx[1] = scale_tSTAT_/init_STAT_;
    dydx[3] = scale_pSTAT_/init_STAT_;
    dydx[4] = scale_tSTAT_/init_STAT_;
    dydx[6] = 2*scale_pSTAT_/init_STAT_;
    dydx[7] = 2*scale_tSTAT_/init_STAT_;
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
