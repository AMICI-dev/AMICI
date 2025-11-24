#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "p.h"
#include "k.h"
#include "xdot.h"

namespace amici {
namespace model_model_steadystate_py {

void xdot_model_steadystate_py(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dx1dt = -2*p1*std::pow(x1, 2) - p2*x1*x2 + 2*p3*x2 + p4*x3 + p5;  // xdot[0]
    dx2dt = p1*std::pow(x1, 2) - p2*x1*x2 - p3*x2 + p4*x3;  // xdot[1]
    dx3dt = -k4*x3 + p2*x1*x2 - p4*x3;  // xdot[2]
}

} // namespace model_model_steadystate_py
} // namespace amici
