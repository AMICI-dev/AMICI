#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_calvetti_py {

void w_model_calvetti_py(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static){
    const realtype V1_ = x[0];
    const realtype V2_ = x[1];
    const realtype V3_ = x[2];
    const realtype f1_ = x[3];
    const realtype f2_ = x[4];
    const realtype f3_ = x[5];
    const realtype V1ss_ = k[0];
    const realtype R1ss_ = k[1];
    const realtype V2ss_ = k[2];
    const realtype R2ss_ = k[3];
    const realtype V3ss_ = k[4];
    const realtype R3ss_ = k[5];
    const realtype Heaviside_0_ = h[0];
    const realtype Heaviside_2_ = h[2];

    realtype &C1ss_ = w[0];
    realtype &L1_ = w[1];
    realtype &L2_ = w[2];
    realtype &L3_ = w[3];
    realtype &p2_ = w[4];
    realtype &p3_ = w[5];
    realtype &s_ = w[6];
    realtype &C2ss_ = w[7];
    realtype &C3ss_ = w[8];
    realtype &R1_ = w[9];
    realtype &R2_ = w[10];
    realtype &R3_ = w[11];
    realtype &f0_ = w[12];
    realtype &rate_of_V1_ = w[13];
    realtype &rate_of_V2_ = w[14];
    realtype &rate_of_V3_ = w[15];
    // static expressions
    if (include_static) {
        C1ss_ = V1ss_/(1.0 - 0.5*R1ss_);  // w[0]
        L1_ = std::pow(R1ss_, 0.33333333333333331)*std::pow(std::fabs(V1ss_), 0.66666666666666663);  // w[1]
        L2_ = std::pow(R2ss_, 0.33333333333333331)*std::pow(std::fabs(V2ss_), 0.66666666666666663);  // w[2]
        L3_ = std::pow(R3ss_, 0.33333333333333331)*std::pow(std::fabs(V3ss_), 0.66666666666666663);  // w[3]
        p2_ = 1.0 - R1ss_;  // w[4]
        p3_ = -R1ss_ - R2ss_ + 1.0;  // w[5]
        C2ss_ = V2ss_/(-0.5*R2ss_ + p2_);  // w[7]
        C3ss_ = V3ss_/(-0.5*R3ss_ + p3_);  // w[8]
    }

    // dynamic expressions
    s_ = Heaviside_0_*Heaviside_2_;  // w[6]
    R1_ = std::pow(L1_, 3)/std::pow(V1_, 2);  // w[9]
    R2_ = std::pow(L2_, 3)/std::pow(V2_, 2);  // w[10]
    R3_ = std::pow(L3_, 3)/std::pow(V3_, 2);  // w[11]
    f0_ = (-R3_*f3_ - f1_*(R1_ + R2_) - f2_*(R2_ + R3_))/R1_ + 2/R1_;  // w[12]
    rate_of_V1_ = -100.0/899.0*V1_/V1ss_ + (1.0/31.0)*s_ + 129.0/899.0 - 2.0/31.0*V1_/(C1ss_*(R3_*f3_ + f1_*(R1_ + R2_) + f2_*(R2_ + R3_)));  // w[13]
    rate_of_V2_ = -100.0/8313.0*V2_/V2ss_ + 151.0/8313.0 - 2.0/163.0*V2_/(C2ss_*(R3_*f3_ + f2_*(R2_ + R3_)));  // w[14]
    rate_of_V3_ = -1.0/121999878.0*V3_/V3ss_ + 500000.0/60999939.0 - 1.0/61.0*V3_/(C3ss_*R3_*f3_);  // w[15]
}

} // namespace model_model_calvetti_py
} // namespace amici
