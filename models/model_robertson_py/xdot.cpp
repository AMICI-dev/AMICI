#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "p.h"
#include "k.h"
#include "dx.h"
#include "xdot.h"

namespace amici {
namespace model_model_robertson_py {

void xdot_model_robertson_py(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const realtype *w){
    de_0 = -dx1dt - p1*x1 + p2*x2*x3;  // xdot[0]
    de_1 = -dx2dt + p1*x1 - p2*x2*x3 - p3*std::pow(x2, 2);  // xdot[1]
    ae_0 = x1 + x2 + x3 - 1;  // xdot[2]
}

} // namespace model_model_robertson_py
} // namespace amici
