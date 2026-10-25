#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_events_py {

void dydx_model_events_py(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    const realtype p4_ = p[3];

    dydx[0] = p4_;
    dydx[1] = p4_;
    dydx[2] = p4_;
}

} // namespace model_model_events_py
} // namespace amici
