#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_nested_events_py {

void dtaudx_model_nested_events_py(realtype *dtaudx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const realtype *tcl, const int ie){
    const realtype Virus_ = x[0];
    const realtype rho_V_ = p[3];
    const realtype delta_V_ = p[4];
    const realtype Heaviside_1_ = h[0];

    switch(ie) {
        case 0:
            dtaudx[0] = 1.0/(Heaviside_1_*Virus_*rho_V_ - Virus_*delta_V_);
            break;
        case 1:
            dtaudx[0] = -1/(-Heaviside_1_*Virus_*rho_V_ + Virus_*delta_V_);
            break;
    }
}

} // namespace model_model_nested_events_py
} // namespace amici
