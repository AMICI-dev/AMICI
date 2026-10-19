#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_dirac_py {

void deltaqB_model_dirac_py(realtype *deltaqB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const int ip, const int ie, const realtype *xdot, const realtype *xdot_eval, const realtype *x_eval, const realtype *xB){
    const realtype dx1dt_ = xdot[0];
    const realtype dx2dt_ = xdot[1];
    const realtype xdot_eval0_ = xdot_eval[0];
    const realtype xdot_eval1_ = xdot_eval[1];
    const realtype xB0_ = xB[0];
    const realtype xB1_ = xB[1];

    switch(ie) {
        case 0:
            switch(ip) {
                case 1:
                    deltaqB[0] = xB0_*(-dx1dt_ + xdot_eval0_) + xB1_*(-dx2dt_ + xdot_eval1_);
                    break;
            }
            break;
    }
}

} // namespace model_model_dirac_py
} // namespace amici
