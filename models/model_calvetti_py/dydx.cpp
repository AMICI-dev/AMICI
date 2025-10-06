#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "k.h"
#include "h.h"
#include "w.h"

namespace amici {
namespace model_model_calvetti_py {

void dydx_model_calvetti_py(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dydx[0] = 1;
    dydx[3] = std::pow(L1, 3)*(2*f1/R1 - 2*(R3*f3 + f1*(R1 + R2) + f2*(R2 + R3))/std::pow(R1, 2) + 4/std::pow(R1, 2))/std::pow(V1, 3);
    dydx[7] = 1;
    dydx[9] = std::pow(L2, 3)*(2*f1 + 2*f2)/(R1*std::pow(V2, 3));
    dydx[14] = 1;
    dydx[15] = std::pow(L3, 3)*(2*f2 + 2*f3)/(R1*std::pow(V3, 3));
    dydx[21] = (-R1 - R2)/R1;
    dydx[22] = 1;
    dydx[27] = (-R2 - R3)/R1;
    dydx[29] = 1;
    dydx[33] = -R3/R1;
}

} // namespace model_model_calvetti_py
} // namespace amici
