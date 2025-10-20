#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"

namespace amici {
namespace model_model_jakstat_adjoint_py {

void dydx_model_jakstat_adjoint_py(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dydx[1] = scale_tSTAT/init_STAT;
    dydx[3] = scale_pSTAT/init_STAT;
    dydx[4] = scale_tSTAT/init_STAT;
    dydx[6] = 2*scale_pSTAT/init_STAT;
    dydx[7] = 2*scale_tSTAT/init_STAT;
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
