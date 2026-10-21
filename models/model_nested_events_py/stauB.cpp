#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_nested_events_py {

void stauB_model_nested_events_py(realtype *stauB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const int ie, const realtype *xdot, const realtype *xdot_eval, const realtype *xdot_pre, const realtype *x_eval, const realtype *xB, const realtype *tcl){
    const realtype dVirusdt_ = xdot[0];
    const realtype xdot_eval0_ = xdot_eval[0];
    const realtype xdot_pre0_ = xdot_pre[0];
    const realtype xB0_ = xB[0];

    switch(ie) {
        case 0:
        case 1:
            stauB[0] = xB0_*(dVirusdt_ - xdot_pre0_);
            break;
        case 2:
            stauB[0] = xB0_*(dVirusdt_ - xdot_eval0_);
            break;
    }
}

} // namespace model_model_nested_events_py
} // namespace amici
