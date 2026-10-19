#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_calvetti_py {

void stauB_model_calvetti_py(realtype *stauB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const int ie, const realtype *xdot, const realtype *xdot_eval, const realtype *xdot_pre, const realtype *x_eval, const realtype *xB, const realtype *tcl){
    const realtype de_0_ = xdot[0];
    const realtype de_1_ = xdot[1];
    const realtype de_2_ = xdot[2];
    const realtype ae_0_ = xdot[3];
    const realtype ae_1_ = xdot[4];
    const realtype ae_2_ = xdot[5];
    const realtype xdot_pre0_ = xdot_pre[0];
    const realtype xdot_pre1_ = xdot_pre[1];
    const realtype xdot_pre2_ = xdot_pre[2];
    const realtype xdot_pre3_ = xdot_pre[3];
    const realtype xdot_pre4_ = xdot_pre[4];
    const realtype xdot_pre5_ = xdot_pre[5];
    const realtype xB0_ = xB[0];
    const realtype xB1_ = xB[1];
    const realtype xB2_ = xB[2];
    const realtype xB3_ = xB[3];
    const realtype xB4_ = xB[4];
    const realtype xB5_ = xB[5];

    switch(ie) {
        case 0:
        case 1:
        case 2:
        case 3:
            stauB[0] = xB0_*(de_0_ - xdot_pre0_) + xB1_*(de_1_ - xdot_pre1_) + xB2_*(de_2_ - xdot_pre2_) + xB3_*(ae_0_ - xdot_pre3_) + xB4_*(ae_1_ - xdot_pre4_) + xB5_*(ae_2_ - xdot_pre5_);
            break;
    }
}

} // namespace model_model_calvetti_py
} // namespace amici
