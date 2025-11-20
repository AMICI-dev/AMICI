#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"
#include "xdot.h"

namespace amici {
namespace model_model_events_py {

void xdot_model_events_py(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dx1dt = -p1*x1*(1 - Heaviside_2);  // xdot[0]
    dx2dt = p2*x1*std::exp(-1.0/10.0*t) - p3*x2;  // xdot[1]
    dx3dt = -Heaviside_4 - x3 + 1;  // xdot[2]
}

} // namespace model_model_events_py
} // namespace amici
