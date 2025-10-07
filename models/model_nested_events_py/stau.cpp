#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "p.h"
#include "h.h"
#include "sx.h"

namespace amici {
namespace model_model_nested_events_py {

void stau_model_nested_events_py(realtype *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const realtype *tcl, const realtype *sx, const int ip, const int ie){
    switch(ie) {
        case 0:
            switch(ip) {
                case 0:
                case 1:
                case 2:
                case 3:
                case 4:
                    stau[0] = sx0/(Heaviside_1*Virus*rho_V - Virus*delta_V);
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
                    stau[0] = -sx0/(-Heaviside_1*Virus*rho_V + Virus*delta_V);
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
