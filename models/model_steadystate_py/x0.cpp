#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_steadystate_py {

void x0_model_steadystate_py(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    const realtype k1_ = k[0];
    const realtype k2_ = k[1];
    const realtype k3_ = k[2];

    x0[0] = k1_;
    x0[1] = k2_;
    x0[2] = k3_;
}

} // namespace model_model_steadystate_py
} // namespace amici
