#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"
#include "xdot.h"
#include "xdot_old.h"
#include "x_old.h"
#include "xB.h"

namespace amici {
namespace model_model_events_py {

void deltaxB_model_events_py(realtype *deltaxB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *x_old, const realtype *xB, const realtype *tcl){
    switch(ie) {
        case 0:
            deltaxB[1] = xB0*(dx1dt - xdot_old0)/(Heaviside_4 + p2*x1*std::exp(-1.0/10.0*t) - p3*x2 + x3 - 1) + xB1*(dx2dt - xdot_old1)/(Heaviside_4 + p2*x1*std::exp(-1.0/10.0*t) - p3*x2 + x3 - 1) + xB2*(dx3dt - xdot_old2)/(Heaviside_4 + p2*x1*std::exp(-1.0/10.0*t) - p3*x2 + x3 - 1);
            deltaxB[2] = -xB0*(dx1dt - xdot_old0)/(Heaviside_4 + p2*x1*std::exp(-1.0/10.0*t) - p3*x2 + x3 - 1) - xB1*(dx2dt - xdot_old1)/(Heaviside_4 + p2*x1*std::exp(-1.0/10.0*t) - p3*x2 + x3 - 1) - xB2*(dx3dt - xdot_old2)/(Heaviside_4 + p2*x1*std::exp(-1.0/10.0*t) - p3*x2 + x3 - 1);
            break;
        case 1:
            deltaxB[0] = xB0*(dx1dt - xdot_old0)/(Heaviside_4 - p1*x1*(1 - Heaviside_2) + x3 - 1) + xB1*(dx2dt - xdot_old1)/(Heaviside_4 - p1*x1*(1 - Heaviside_2) + x3 - 1) + xB2*(dx3dt - xdot_old2)/(Heaviside_4 - p1*x1*(1 - Heaviside_2) + x3 - 1);
            deltaxB[2] = -xB0*(dx1dt - xdot_old0)/(Heaviside_4 - p1*x1*(1 - Heaviside_2) + x3 - 1) - xB1*(dx2dt - xdot_old1)/(Heaviside_4 - p1*x1*(1 - Heaviside_2) + x3 - 1) - xB2*(dx3dt - xdot_old2)/(Heaviside_4 - p1*x1*(1 - Heaviside_2) + x3 - 1);
            break;
    }
}

} // namespace model_model_events_py
} // namespace amici
