#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_nested_events_py {

void xdot_model_nested_events_py(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    const realtype Virus_ = x[0];
    const realtype rho_V_ = p[3];
    const realtype delta_V_ = p[4];
    const realtype Heaviside_1_ = h[0];

    realtype &dVirusdt_ = xdot[0];
    dVirusdt_ = Heaviside_1_*Virus_*rho_V_ - Virus_*delta_V_;  // xdot[0]
}

} // namespace model_model_nested_events_py
} // namespace amici
