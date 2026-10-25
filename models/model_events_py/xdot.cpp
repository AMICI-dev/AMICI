#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_events_py {

void xdot_model_events_py(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    const realtype x1_ = x[0];
    const realtype x2_ = x[1];
    const realtype x3_ = x[2];
    const realtype p1_ = p[0];
    const realtype p2_ = p[1];
    const realtype p3_ = p[2];
    const realtype Heaviside_2_ = h[2];
    const realtype Heaviside_4_ = h[4];

    realtype &dx1dt_ = xdot[0];
    realtype &dx2dt_ = xdot[1];
    realtype &dx3dt_ = xdot[2];
    dx1dt_ = -p1_*x1_*(1 - Heaviside_2_);  // xdot[0]
    dx2dt_ = p2_*x1_*std::exp(-1.0/10.0*t) - p3_*x2_;  // xdot[1]
    dx3dt_ = -Heaviside_4_ - x3_ + 1;  // xdot[2]
}

} // namespace model_model_events_py
} // namespace amici
