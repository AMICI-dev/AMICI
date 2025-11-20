#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "p.h"
#include "h.h"
#include "xdot.h"

namespace amici {
namespace model_model_nested_events_py {

void xdot_model_nested_events_py(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dVirusdt = Heaviside_1*Virus*rho_V - Virus*delta_V;  // xdot[0]
}

} // namespace model_model_nested_events_py
} // namespace amici
