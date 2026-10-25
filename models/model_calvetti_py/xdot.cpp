#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_calvetti_py {

void xdot_model_calvetti_py(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const realtype *w){
    const realtype f1_ = x[3];
    const realtype f2_ = x[4];
    const realtype f3_ = x[5];
    const realtype dV1dt_ = dx[0];
    const realtype dV2dt_ = dx[1];
    const realtype dV3dt_ = dx[2];
    const realtype f0_ = w[12];
    const realtype rate_of_V1_ = w[13];
    const realtype rate_of_V2_ = w[14];
    const realtype rate_of_V3_ = w[15];

    realtype &de_0_ = xdot[0];
    realtype &de_1_ = xdot[1];
    realtype &de_2_ = xdot[2];
    realtype &ae_0_ = xdot[3];
    realtype &ae_1_ = xdot[4];
    realtype &ae_2_ = xdot[5];
    de_0_ = -dV1dt_ + rate_of_V1_;  // xdot[0]
    de_1_ = -dV2dt_ + rate_of_V2_;  // xdot[1]
    de_2_ = -dV3dt_ + rate_of_V3_;  // xdot[2]
    ae_0_ = f0_ - f1_ - rate_of_V1_;  // xdot[3]
    ae_1_ = f1_ - f2_ - rate_of_V2_;  // xdot[4]
    ae_2_ = f2_ - f3_ - rate_of_V3_;  // xdot[5]
}

} // namespace model_model_calvetti_py
} // namespace amici
