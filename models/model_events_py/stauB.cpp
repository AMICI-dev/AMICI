#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_events_py {

void stauB_model_events_py(realtype *stauB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const int ie, const realtype *xdot, const realtype *xdot_eval, const realtype *xdot_pre, const realtype *x_eval, const realtype *xB, const realtype *tcl){
    const realtype dx1dt_ = xdot[0];
    const realtype dx2dt_ = xdot[1];
    const realtype dx3dt_ = xdot[2];
    const realtype xdot_pre0_ = xdot_pre[0];
    const realtype xdot_pre1_ = xdot_pre[1];
    const realtype xdot_pre2_ = xdot_pre[2];
    const realtype xB0_ = xB[0];
    const realtype xB1_ = xB[1];
    const realtype xB2_ = xB[2];

    switch(ie) {
        case 0:
        case 1:
        case 2:
        case 3:
        case 4:
        case 5:
            stauB[0] = xB0_*(dx1dt_ - xdot_pre0_) + xB1_*(dx2dt_ - xdot_pre1_) + xB2_*(dx3dt_ - xdot_pre2_);
            break;
    }
}

} // namespace model_model_events_py
} // namespace amici
