#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_neuron_py {

void stauB_model_neuron_py(realtype *stauB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const int ie, const realtype *xdot, const realtype *xdot_eval, const realtype *xdot_pre, const realtype *x_eval, const realtype *xB, const realtype *tcl){
    const realtype dvdt_ = xdot[0];
    const realtype dudt_ = xdot[1];
    const realtype xdot_eval1_ = xdot_eval[1];
    const realtype xB0_ = xB[0];
    const realtype xB1_ = xB[1];

    switch(ie) {
        case 0:
            stauB[0] = dvdt_*xB0_ + xB1_*(dudt_ - xdot_eval1_);
            break;
    }
}

} // namespace model_model_neuron_py
} // namespace amici
