#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_calvetti_py {

void dydx_model_calvetti_py(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    const realtype V1_ = x[0];
    const realtype V2_ = x[1];
    const realtype V3_ = x[2];
    const realtype f1_ = x[3];
    const realtype f2_ = x[4];
    const realtype f3_ = x[5];
    const realtype L1_ = w[1];
    const realtype L2_ = w[2];
    const realtype L3_ = w[3];
    const realtype R1_ = w[9];
    const realtype R2_ = w[10];
    const realtype R3_ = w[11];

    dydx[0] = 1;
    dydx[3] = std::pow(L1_, 3)*(2*f1_/R1_ - 2*(R3_*f3_ + f1_*(R1_ + R2_) + f2_*(R2_ + R3_))/std::pow(R1_, 2) + 4/std::pow(R1_, 2))/std::pow(V1_, 3);
    dydx[7] = 1;
    dydx[9] = std::pow(L2_, 3)*(2*f1_ + 2*f2_)/(R1_*std::pow(V2_, 3));
    dydx[14] = 1;
    dydx[15] = std::pow(L3_, 3)*(2*f2_ + 2*f3_)/(R1_*std::pow(V3_, 3));
    dydx[21] = (-R1_ - R2_)/R1_;
    dydx[22] = 1;
    dydx[27] = (-R2_ - R3_)/R1_;
    dydx[29] = 1;
    dydx[33] = -R3_/R1_;
}

} // namespace model_model_calvetti_py
} // namespace amici
