#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_nested_events_py {

void deltaxB_model_nested_events_py(realtype *deltaxB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const int ie, const realtype *xdot, const realtype *xdot_eval, const realtype *x_eval, const realtype *xB, const realtype *tcl){
    const realtype rho_V_ = p[3];
    const realtype delta_V_ = p[4];
    const realtype Heaviside_1_ = h[0];
    const realtype dVirusdt_ = xdot[0];
    const realtype xdot_eval0_ = xdot_eval[0];
    const realtype x_eval0_ = x_eval[0];
    const realtype xB0_ = xB[0];

    switch(ie) {
        case 0:
            deltaxB[0] = xB0_*(dVirusdt_ - xdot_eval0_)/(Heaviside_1_*rho_V_*x_eval0_ - delta_V_*x_eval0_);
            break;
        case 1:
            deltaxB[0] = -xB0_*(dVirusdt_ - xdot_eval0_)/(-Heaviside_1_*rho_V_*x_eval0_ + delta_V_*x_eval0_);
            break;
    }
}

} // namespace model_model_nested_events_py
} // namespace amici
