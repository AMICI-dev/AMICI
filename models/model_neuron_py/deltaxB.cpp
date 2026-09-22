#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_neuron_py {

void deltaxB_model_neuron_py(realtype *deltaxB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *x_old, const realtype *xB, const realtype *tcl){
    const realtype I0_ = k[1];
    const realtype dvdt_ = xdot[0];
    const realtype dudt_ = xdot[1];
    const realtype xdot_old0_ = xdot_old[0];
    const realtype xdot_old1_ = xdot_old[1];
    const realtype x_old0_ = x_old[0];
    const realtype x_old1_ = x_old[1];
    const realtype xB0_ = xB[0];
    const realtype xB1_ = xB[1];

    switch(ie) {
        case 0:
            deltaxB[0] = xB0_*(xdot_old0_/(I0_ + (1.0/25.0)*std::pow(x_old0_, 2) + 5*x_old0_ - x_old1_ + 140) + (dvdt_ - xdot_old0_)/(I0_ + (1.0/25.0)*std::pow(x_old0_, 2) + 5*x_old0_ - x_old1_ + 140) - 1) + xB1_*(dudt_ - xdot_old1_)/(I0_ + (1.0/25.0)*std::pow(x_old0_, 2) + 5*x_old0_ - x_old1_ + 140);
            break;
    }
}

} // namespace model_model_neuron_py
} // namespace amici
