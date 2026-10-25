#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_nested_events_py {

void stau_model_nested_events_py(realtype *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const realtype *tcl, const realtype *sx, const int ip, const int ie){
    const realtype Virus_ = x[0];
    const realtype rho_V_ = p[3];
    const realtype delta_V_ = p[4];
    const realtype Heaviside_1_ = h[0];
    const realtype sx0_ = sx[0];

    switch(ie) {
        case 0:
            switch(ip) {
                case 0:
                case 1:
                case 2:
                case 3:
                case 4:
                    stau[0] = sx0_/(Heaviside_1_*Virus_*rho_V_ - Virus_*delta_V_);
                    break;
            }
            break;
        case 1:
            switch(ip) {
                case 0:
                case 1:
                case 2:
                case 3:
                case 4:
                    stau[0] = -sx0_/(-Heaviside_1_*Virus_*rho_V_ + Virus_*delta_V_);
                    break;
            }
            break;
        case 2:
            switch(ip) {
                case 2:
                    stau[0] = -1;
                    break;
            }
            break;
    }
}

} // namespace model_model_nested_events_py
} // namespace amici
