#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "p.h"
#include "h.h"
#include "xdot.h"

namespace amici {
namespace model_model_dirac_py {

void xdot_model_dirac_py(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dx1dt = -p1*x1;  // xdot[0]
    dx2dt = p3*x1 - p4*x2;  // xdot[1]
}

} // namespace model_model_dirac_py
} // namespace amici
