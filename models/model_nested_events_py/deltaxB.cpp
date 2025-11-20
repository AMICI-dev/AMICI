#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "p.h"
#include "h.h"
#include "xdot.h"
#include "xdot_old.h"
#include "x_old.h"
#include "xB.h"

namespace amici {
namespace model_model_nested_events_py {

void deltaxB_model_nested_events_py(realtype *deltaxB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *x_old, const realtype *xB, const realtype *tcl){
    switch(ie) {
        case 0:
            deltaxB[0] = xB0*(dVirusdt - xdot_old0)/(Heaviside_1*Virus*rho_V - Virus*delta_V);
            break;
        case 1:
            deltaxB[0] = -xB0*(dVirusdt - xdot_old0)/(-Heaviside_1*Virus*rho_V + Virus*delta_V);
            break;
    }
}

} // namespace model_model_nested_events_py
} // namespace amici
