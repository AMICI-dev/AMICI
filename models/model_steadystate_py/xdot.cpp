#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_steadystate_py {

void xdot_model_steadystate_py(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    const realtype x1_ = x[0];
    const realtype x2_ = x[1];
    const realtype x3_ = x[2];
    const realtype p1_ = p[0];
    const realtype p2_ = p[1];
    const realtype p3_ = p[2];
    const realtype p4_ = p[3];
    const realtype p5_ = p[4];
    const realtype k4_ = k[3];

    realtype &dx1dt_ = xdot[0];
    realtype &dx2dt_ = xdot[1];
    realtype &dx3dt_ = xdot[2];
    dx1dt_ = -2*p1_*std::pow(x1_, 2) - p2_*x1_*x2_ + 2*p3_*x2_ + p4_*x3_ + p5_;  // xdot[0]
    dx2dt_ = p1_*std::pow(x1_, 2) - p2_*x1_*x2_ - p3_*x2_ + p4_*x3_;  // xdot[1]
    dx3dt_ = -k4_*x3_ + p2_*x1_*x2_ - p4_*x3_;  // xdot[2]
}

} // namespace model_model_steadystate_py
} // namespace amici
