#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_nested_events_py {

void deltaxB_model_nested_events_py(realtype *deltaxB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *x_old, const realtype *xB, const realtype *tcl){
    const realtype rho_V_ = p[3];
    const realtype delta_V_ = p[4];
    const realtype Heaviside_1_ = h[0];
    const realtype dVirusdt_ = xdot[0];
    const realtype xdot_old0_ = xdot_old[0];
    const realtype x_old0_ = x_old[0];
    const realtype xB0_ = xB[0];

    switch(ie) {
        case 0:
            deltaxB[0] = xB0_*(dVirusdt_ - xdot_old0_)/(Heaviside_1_*rho_V_*x_old0_ - delta_V_*x_old0_);
            break;
        case 1:
            deltaxB[0] = -xB0_*(dVirusdt_ - xdot_old0_)/(-Heaviside_1_*rho_V_*x_old0_ + delta_V_*x_old0_);
            break;
    }
}

} // namespace model_model_nested_events_py
} // namespace amici
