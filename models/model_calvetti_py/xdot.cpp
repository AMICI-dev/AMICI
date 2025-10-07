#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "k.h"
#include "h.h"
#include "dx.h"
#include "w.h"
#include "xdot.h"

namespace amici {
namespace model_model_calvetti_py {

void xdot_model_calvetti_py(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const realtype *w){
    de_0 = -dV1dt + rate_of_V1;  // xdot[0]
    de_1 = -dV2dt + rate_of_V2;  // xdot[1]
    de_2 = -dV3dt + rate_of_V3;  // xdot[2]
    ae_0 = f0 - f1 - rate_of_V1;  // xdot[3]
    ae_1 = f1 - f2 - rate_of_V2;  // xdot[4]
    ae_2 = f2 - f3 - rate_of_V3;  // xdot[5]
}

} // namespace model_model_calvetti_py
} // namespace amici
