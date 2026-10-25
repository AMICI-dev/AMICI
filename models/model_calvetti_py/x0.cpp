#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_calvetti_py {

void x0_model_calvetti_py(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    const realtype V1ss_ = k[0];
    const realtype V2ss_ = k[2];
    const realtype V3ss_ = k[4];

    x0[0] = V1ss_;
    x0[1] = V2ss_;
    x0[2] = V3ss_;
    x0[3] = 1.0;
    x0[4] = 1.0;
    x0[5] = 1.0;
}

} // namespace model_model_calvetti_py
} // namespace amici
