#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "k.h"
#include "h.h"
#include "w.h"

namespace amici {
namespace model_model_calvetti_py {

void w_model_calvetti_py(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        C1ss = V1ss/(1.0 - 0.5*R1ss);  // w[0]
        C2ss = V2ss/(-R1ss - 0.5*R2ss + 1.0);  // w[1]
        C3ss = V3ss/(-R1ss - R2ss - 0.5*R3ss + 1.0);  // w[2]
        L1 = std::pow(R1ss, 0.33333333333333331)*std::pow(std::fabs(V1ss), 0.66666666666666663);  // w[3]
        L2 = std::pow(R2ss, 0.33333333333333331)*std::pow(std::fabs(V2ss), 0.66666666666666663);  // w[4]
        L3 = std::pow(R3ss, 0.33333333333333331)*std::pow(std::fabs(V3ss), 0.66666666666666663);  // w[5]
        p2 = 1.0 - R1ss;  // w[6]
        p3 = -R1ss - R2ss + 1.0;  // w[7]
    }

    // dynamic expressions
    s = Heaviside_0*Heaviside_2;  // w[8]
    R1 = std::pow(L1, 3)/std::pow(V1, 2);  // w[9]
    R2 = std::pow(L2, 3)/std::pow(V2, 2);  // w[10]
    R3 = std::pow(L3, 3)/std::pow(V3, 2);  // w[11]
    f0 = (-R3*f3 - f1*(R1 + R2) - f2*(R2 + R3))/R1 + 2/R1;  // w[12]
    rate_of_V1 = -100.0/899.0*V1/V1ss + (1.0/31.0)*s + 129.0/899.0 - 2.0/31.0*V1/(C1ss*(R3*f3 + f1*(R1 + R2) + f2*(R2 + R3)));  // w[13]
    rate_of_V2 = -100.0/8313.0*V2/V2ss + 151.0/8313.0 - 2.0/163.0*V2/(C2ss*(R3*f3 + f2*(R2 + R3)));  // w[14]
    rate_of_V3 = -1.0/121999878.0*V3/V3ss + 500000.0/60999939.0 - 1.0/61.0*V3/(C3ss*R3*f3);  // w[15]
}

} // namespace model_model_calvetti_py
} // namespace amici
