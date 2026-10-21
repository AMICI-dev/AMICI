#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_neuron_py {

void deltasx_model_neuron_py(realtype *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const int ip, const int ie, const realtype *xdot, const realtype *xdot_eval, const realtype *xdot_pre, const realtype *sx, const realtype *sx_eval, const realtype *stau, const realtype *tcl, const realtype *x_eval){
    const realtype dvdt_ = xdot[0];
    const realtype dudt_ = xdot[1];
    const realtype xdot_eval1_ = xdot_eval[1];
    const realtype xdot_pre0_ = xdot_pre[0];
    const realtype xdot_pre1_ = xdot_pre[1];
    const realtype sx0_ = sx[0];
    const realtype sx1_ = sx[1];
    const realtype sx_eval1_ = sx_eval[1];
    const realtype stau0_ = stau[0];

    switch(ie) {
        case 0:
            switch(ip) {
                case 0:
                case 1:
                    deltasx[0] = stau0_*xdot_pre0_ + stau0_*(dvdt_ - xdot_pre0_) - sx0_;
                    deltasx[1] = -stau0_*xdot_eval1_ + stau0_*xdot_pre1_ + stau0_*(dudt_ - xdot_pre1_) - sx1_ + sx_eval1_;
                    break;
                case 2:
                    deltasx[0] = stau0_*xdot_pre0_ + stau0_*(dvdt_ - xdot_pre0_) - sx0_ - 1;
                    deltasx[1] = -stau0_*xdot_eval1_ + stau0_*xdot_pre1_ + stau0_*(dudt_ - xdot_pre1_) - sx1_ + sx_eval1_;
                    break;
                case 3:
                    deltasx[0] = stau0_*xdot_pre0_ + stau0_*(dvdt_ - xdot_pre0_) - sx0_;
                    deltasx[1] = -stau0_*xdot_eval1_ + stau0_*xdot_pre1_ + stau0_*(dudt_ - xdot_pre1_) - sx1_ + sx_eval1_ + 1;
                    break;
            }
            break;
    }
}

} // namespace model_model_neuron_py
} // namespace amici
