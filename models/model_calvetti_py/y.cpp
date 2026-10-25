#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_calvetti_py {

void y_model_calvetti_py(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    const realtype V1_ = x[0];
    const realtype V2_ = x[1];
    const realtype V3_ = x[2];
    const realtype f1_ = x[3];
    const realtype f2_ = x[4];
    const realtype f0_ = w[12];

    y[0] = V1_;
    y[1] = V2_;
    y[2] = V3_;
    y[3] = f0_;
    y[4] = f1_;
    y[5] = f2_;
}

} // namespace model_model_calvetti_py
} // namespace amici
