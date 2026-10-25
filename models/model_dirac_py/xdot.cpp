#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_dirac_py {

void xdot_model_dirac_py(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    const realtype x1_ = x[0];
    const realtype x2_ = x[1];
    const realtype p1_ = p[0];
    const realtype p3_ = p[2];
    const realtype p4_ = p[3];

    realtype &dx1dt_ = xdot[0];
    realtype &dx2dt_ = xdot[1];
    dx1dt_ = -p1_*x1_;  // xdot[0]
    dx2dt_ = p3_*x1_ - p4_*x2_;  // xdot[1]
}

} // namespace model_model_dirac_py
} // namespace amici
