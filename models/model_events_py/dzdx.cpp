#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_events_py {

void dzdx_model_events_py(realtype *dzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h){
    const realtype x1_ = x[0];
    const realtype x2_ = x[1];
    const realtype x3_ = x[2];
    const realtype p1_ = p[0];
    const realtype p2_ = p[1];
    const realtype p3_ = p[2];
    const realtype Heaviside_2_ = h[2];
    const realtype Heaviside_4_ = h[4];

    switch(ie) {
        case 0:
            dzdx[2] = -1/(Heaviside_4_ + p2_*x1_*std::exp(-1.0/10.0*t) - p3_*x2_ + x3_ - 1);
            dzdx[4] = 1.0/(Heaviside_4_ + p2_*x1_*std::exp(-1.0/10.0*t) - p3_*x2_ + x3_ - 1);
            break;
        case 1:
            dzdx[1] = -1/(Heaviside_4_ - p1_*x1_*(1 - Heaviside_2_) + x3_ - 1);
            dzdx[5] = 1.0/(Heaviside_4_ - p1_*x1_*(1 - Heaviside_2_) + x3_ - 1);
            break;
    }
}

} // namespace model_model_events_py
} // namespace amici
