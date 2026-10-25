#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_robertson_py {

void xdot_model_robertson_py(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const realtype *w){
    const realtype x1_ = x[0];
    const realtype x2_ = x[1];
    const realtype x3_ = x[2];
    const realtype p1_ = p[0];
    const realtype p2_ = p[1];
    const realtype p3_ = p[2];
    const realtype dx1dt_ = dx[0];
    const realtype dx2dt_ = dx[1];

    realtype &de_0_ = xdot[0];
    realtype &de_1_ = xdot[1];
    realtype &ae_0_ = xdot[2];
    de_0_ = -dx1dt_ - p1_*x1_ + p2_*x2_*x3_;  // xdot[0]
    de_1_ = -dx2dt_ + p1_*x1_ - p2_*x2_*x3_ - p3_*std::pow(x2_, 2);  // xdot[1]
    ae_0_ = x1_ + x2_ + x3_ - 1;  // xdot[2]
}

} // namespace model_model_robertson_py
} // namespace amici
