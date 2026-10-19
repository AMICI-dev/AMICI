#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_neuron_py {

void deltaxB_model_neuron_py(realtype *deltaxB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const int ie, const realtype *xdot, const realtype *xdot_eval, const realtype *x_eval, const realtype *xB, const realtype *tcl){
    const realtype I0_ = k[1];
    const realtype dvdt_ = xdot[0];
    const realtype dudt_ = xdot[1];
    const realtype xdot_eval0_ = xdot_eval[0];
    const realtype xdot_eval1_ = xdot_eval[1];
    const realtype x_eval0_ = x_eval[0];
    const realtype x_eval1_ = x_eval[1];
    const realtype xB0_ = xB[0];
    const realtype xB1_ = xB[1];

    switch(ie) {
        case 0:
            deltaxB[0] = xB0_*(xdot_eval0_/(I0_ + (1.0/25.0)*std::pow(x_eval0_, 2) + 5*x_eval0_ - x_eval1_ + 140) + (dvdt_ - xdot_eval0_)/(I0_ + (1.0/25.0)*std::pow(x_eval0_, 2) + 5*x_eval0_ - x_eval1_ + 140) - 1) + xB1_*(dudt_ - xdot_eval1_)/(I0_ + (1.0/25.0)*std::pow(x_eval0_, 2) + 5*x_eval0_ - x_eval1_ + 140);
            break;
    }
}

} // namespace model_model_neuron_py
} // namespace amici
