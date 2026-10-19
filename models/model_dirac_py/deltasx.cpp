#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_dirac_py {

void deltasx_model_dirac_py(realtype *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const int ip, const int ie, const realtype *xdot, const realtype *xdot_eval, const realtype *sx, const realtype *stau, const realtype *tcl, const realtype *x_eval){
    const realtype dx1dt_ = xdot[0];
    const realtype dx2dt_ = xdot[1];
    const realtype xdot_eval0_ = xdot_eval[0];
    const realtype xdot_eval1_ = xdot_eval[1];
    const realtype stau0_ = stau[0];

    switch(ie) {
        case 0:
            switch(ip) {
                case 0:
                case 1:
                case 2:
                case 3:
                    deltasx[0] = stau0_*(dx1dt_ - xdot_eval0_);
                    deltasx[1] = stau0_*(dx2dt_ - xdot_eval1_);
                    break;
            }
            break;
    }
}

} // namespace model_model_dirac_py
} // namespace amici
